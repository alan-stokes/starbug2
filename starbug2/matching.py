"""
Starbug matching functions
Primarily this is the main routines for dither/band/generic matching which are
 at the core of starbug2 and starbug2-match
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, hstack, Column, vstack

import starbug2
from starbug2.constants import VERBOSE_TAG, CAT_NUM
from starbug2.param import load_params
from starbug2.utils import (
    Loading, printf, remove_duplicates, p_error, fill_nan, tab2array,
    find_col_names, flux2mag, h_cascade, warn, puts)

# keys for catalogue fields.
_OBS = "OBSERVTN"
_FILTER = "FILTER"
_VISIT = "VISIT"
_DETECTOR = "DETECTOR"
_EXPOSURE = "EXPOSURE"

class GenericMatch(object):

    @staticmethod
    def mask_catalogues(catalogues, mask):
        """ takes catalogues and masks and removes catalogues which don't
        match the mask.

        :param catalogues: the catalogues to match.
        :type catalogues: list (astropy.Tables)
        :param mask: the mask to apply.
        :return: an astro table with masked catalogues.
        :rtype astropy.Table
        """
        masked = Table(None)

        if mask is None or type(mask) not in (list,np.ndarray):
            return masked
        if len(catalogues) != len(catalogues):
            return masked

        for subset, cat in zip(mask,catalogues):
            if subset is not None:
                if type(subset) == list:
                    subset=np.array(subset)
                if len(subset) == len(cat):
                    masked=vstack( (masked,cat[~subset]) )
                    cat.remove_rows(~subset)
        return masked


    def __init__(
            self, threshold=None, col_names=None, filter_string=None,
            verbose=None, p_file=None, method="Generic Matching",
            load=Loading(1)):
        """
        constructor for the generic match.

        :param threshold: Separation threshold in arc-seconds
        :type threshold: float or None
        :param col_names: List of str column names to include in the matching.
        Everything else will be discarded.
        :type col_names: list or None
        :param filter_string: Specifically set the filter of the catalogues
        :type filter_string: str or None
        :param verbose: Include verbose outputs
        :type verbose: int or None
        :param p_file: Parameter filename
        :type p_file: str or None
        :param load: the loading object
        :type load: Loading
        """
        options=load_params(p_file)
        self.threshold = options.get("MATCH_THRESH")
        self.filter = options.get(_FILTER)
        self.verbose = options.get(VERBOSE_TAG)
        self.method = method

        if threshold is not None:
            self.threshold = threshold
        self.threshold *= u.arcsec

        if filter_string is not None:
            self.filter = filter_string
        if verbose is not None:
            self.verbose = verbose

        self.col_names = col_names
        self.load = load
        
    def log(self,msg):
        """
        logs messages only when in verbose mode.
        :param msg: message to log
        :return: None
        """
        if self.verbose:
            printf(msg)

    def __str__(self):
        """
        string representation fo the generic match class.
        :return: str
        """
        s=[ "%s:" % self.method,
            "Filter: %s" % self.filter,
            "Col names: %s" % self.col_names,
            "Threshold: %s\"" % self.threshold]
        return "\n".join(s)

    def __call__(self, *args, **kwargs):
        """
        main entrance method.

        :param args: args to the matching algorithm.
        :param kwargs: extra args.
        :return: matched catalogue
        :rtype: astropy.Table
        """
        return self.match(*args, **kwargs)

    def init_catalogues(self, catalogues):
        """
        This function is a bit of a "do everything" function

        It takes the input catalogues and removes any columns that aren't
        included in the self.colnames list. If this is None, then all columns
        are kept. Additionally, it initialises the loading bar with the summed
        length of all the input catalogues.
        Finally, it attempts to set the photometric filter that is being used

        :param catalogues: The input catalogues to work on
        :type catalogues: list (astropy.Tables)
        :return: The cleaned list of input catalogues
        :rtype: list (astropy.Table)
        """

        ## Must copy here maybe?
        if len(catalogues)>=2:
            self.load=Loading(
                sum(len(cat) for cat in catalogues[1:]), msg="initialising")
            if self.verbose: self.load.show()

        # initialise the column names if it wasn't already set
        if self.col_names is None:
            self.col_names = []
            for cat in catalogues:
                self.col_names += cat.col_names
        self.col_names = remove_duplicates(self.col_names)

        # clean out the column names not included in self.col_names
        for n,catalogue in enumerate(catalogues):
            keep= set(catalogue.col_names) & set(self.col_names)
            keep=sorted(keep, key= lambda s:self.col_names.index(s))
            catalogues[n]=catalogue[keep]

        # Attempt to get a value for filter if not already set
        if not self.filter:
            if (filter_string := catalogues[0].meta.get(_FILTER)) is None:
                filter_string = "MAG"
            self.filter = filter_string

        return catalogues


    def match(self, catalogues, join_type="or", mask=None, cartesian=False, **kwargs):
        """
        This matching works as a basic match. Everything is included and the column
        names have _N appended to the end. 

        Parameters
        ----------
        join_type : str
            Joing method
            "or" include sources in any catalogue
            "and" only include sources in all catalogues

        mask : list
            in prep. 
        
        Returns
        -------
        base : astropy.Table
            Matched catalogue.
        """
        catalogues=self.init_catalogues(catalogues)
        if CAT_NUM in self.col_names: self.col_names.remove(CAT_NUM)
        masked= self.mask_catalogues(catalogues, mask)
        base=self.build_meta(catalogues)

        if join_type=="and": p_error("join_type 'and' not fully implemented\n")

        for n,cat in enumerate(catalogues,1): # Bulk matching processes (column naming)
            self.load.msg="matching: %d"%n
            tmp=self._match(base,cat, join_type=join_type, cartesian=cartesian)
            tmp.rename_columns( tmp.colnames, ["%s_%d"%(name,n) for name in tmp.colnames] )
            base=fill_nan(hstack((base,tmp)))

        if len(masked): # Add in any masked bits
            masked.rename_columns( masked.colnames, ["%s_0"%n for n in masked.colnames])
            base=fill_nan(vstack(( base,masked)))
        return base

    def _match(self, base, cat, join_type="or", cartesian=False):
        """
        Base matching function between two catalogues
        
        Parameters
        ----------
        cat1 and cat2 : `astropy.table.Table`
            two astropy tables containing columns with "RA/DEC" in the column names.
            This could be RA_1, RA_2 .. 
            If several columns are located, they will be nanmeaned together

        Returns
        -------
        idx,d2d,d3d : 
            the same as SkyCoord.match_to_catalog_3d
        """
        if not len(base): return cat.copy()

        base=fill_nan(base.copy())
        colnames=[n for n in self.col_names if n in cat.col_names]
        cat=fill_nan(cat[colnames].copy())

        if not cartesian:
            _ra_cols= list(name for name in base.col_names if "RA" in name)
            _dec_cols= list(name for name in base.col_names if "DEC" in name)
            _ra= np.nanmean(tab2array(base, col_names=_ra_cols), axis=1)
            _dec=np.nanmean(tab2array(base, col_names=_dec_cols), axis=1)
            skycoord1=SkyCoord( ra=_ra*u.deg, dec=_dec*u.deg)

            _ra_cols= list(name for name in cat.col_names if "RA" in name)
            _dec_cols= list(name for name in cat.col_names if "DEC" in name)
            _ra= np.nanmean(tab2array(cat, col_names=_ra_cols), axis=1)
            _dec=np.nanmean(tab2array(cat, col_names=_dec_cols), axis=1)
            skycoord2=SkyCoord( ra=_ra*u.deg, dec=_dec*u.deg)
        else:
            _x_cols= list(name for name in base.col_names if name[0] == "x")
            _y_cols= list(name for name in base.col_names if name[0] == "y")
            _x= np.nanmean(tab2array(base, col_names=_x_cols), axis=1)
            _y=np.nanmean(tab2array(base, col_names=_y_cols), axis=1)
            skycoord1=SkyCoord( x=_x, y=_y, z=np.zeros(len(_x)), representation_type="cartesian")

            _x_cols= list(name for name in cat.col_names if name[0] == "x")
            _y_cols= list(name for name in cat.col_names if name[0] == "y")
            _x= np.nanmean(tab2array(cat, col_names=_x_cols), axis=1)
            _y=np.nanmean(tab2array(cat, col_names=_y_cols), axis=1)
            skycoord2=SkyCoord( x=_x, y=_y, z=np.zeros(len(_x)), representation_type="cartesian")


        #######################
        # The actual Matching #
        #######################
        idx,d2d,d3d=skycoord2.match_to_catalog_3d(skycoord1)
        tmp=Table(np.full( (len(base),len(colnames)),np.nan), names=colnames, dtype=cat[colnames].dtype)

        if cartesian:
            dist=d3d
            threshold=self.threshold.value
        else: 
            dist=d2d
            threshold=self.threshold

        for src,IDX,sep in zip(cat, idx, dist):
            self.load()
            if self.verbose: self.load.show()

            if (sep<=threshold) and (sep==min(dist[idx==IDX])): ##GOODMATCH
                tmp[IDX]=src
            elif join_type=="or": ## Append a source
                tmp.add_row(src)
        
        return tmp


    def finish_matching(self, tab, error_column="eflux", num_thresh=-1, discard_outliers=False, zpmag=0, colnames=None):
        """
        Averaging all the values. Combining source flags and building a NUM column

        Parameters
        ----------
        tab : `astropy.table.Table`
            Table to work on

        error_column : str
            Column containing resultant photometric errors to be used to calculate the magnitude error
            "eflux" - use the eflux column (the normal photometric error column)
            "stdflux"   - use "stdflux" column as error on flux

        num_thresh : int
            Minimum number of matches a source must have.
            NUM values smaller than this will be removed from the table.
            If num_thresh<=0, no cropping will happen

        discard_outliers : bool
            Choose whether to remove outling values from the averaging. 
            This may be usful if some flux values are wildly different from others in the set

        zpmag : float
            Zero point (Magnitude) to be applied to the magnitude after it is calculated

        colnames : list
            List of colnames to include in the averaging. If None, use self.colnames instead

        Returns
        -------
        av : astropy.Table
            An averaged version of the input table
        """
        flags=np.full(len(tab),starbug2.SRC_GOOD, dtype=np.uint16)
        av=Table(None)#np.full((len(tab),len(self.colnames)),np.nan), names=self.colnames)
        
        if colnames is None: colnames=self.col_names
        for ii,name in enumerate(colnames):
            #print(name,av.colnames)
            if (all_cols:=find_col_names(tab, name)):
                col=Column(None, name=name)
                ar=tab2array(tab, col_names=all_cols)
                if ar.shape[1]>1:
                    if name=="flux":
                        col=Column(np.nanmedian(ar,axis=1), name=name)
                        mean=np.nanmean(ar,axis=1)
                        if "stdflux" not in self.col_names:
                            av.add_column(Column(np.nanstd(ar,axis=1),name="stdflux"),index=ii+1)
                        ## if median and mean are >5% different, flag as SRC_VAR
                        flags[ np.abs(mean-col)>(col/5.0)] |= starbug2.SRC_VAR
                    elif name== "eflux":
                        col=Column(np.sqrt(np.nansum(ar*ar, axis=1)), name=name)
                    elif name=="stdflux": 
                        col=Column(np.nanmedian(ar,axis=1),name=name)
                    elif name=="flag":
                        col=Column(flags, name=name)
                        for fcol in ar.T: flags|=fcol.astype(np.uint16)
                    elif name=="NUM":
                        col=Column(np.nansum(ar, axis=1), name=name)
                    elif name==CAT_NUM:
                        col=Column(all_cols[0],name=name)
                    else:
                        col=Column(np.nanmedian(ar, axis=1),name=name)
                else:
                    col=tab[all_cols[0]]
                    col.name=name
                
                #av[name]=col
                av.add_column(col,index=ii)
            #else: av.remove_column(name) ## Clean empty columns in table

        av["flag"]=Column(flags,name="flag")
        if "flux" in av.colnames:
            ecol=av[error_column] if error_column in av.colnames else None
            mag,magerr=flux2mag(av["flux"], flux_err=ecol)
            mag+=zpmag

            if self.filter in av.colnames: av.remove_column(self.filter)
            if "e%s"%self.filter in av.colnames: av.remove_column("e%s"%self.filter)
            av.add_column(mag,name=self.filter)
            av.add_column(magerr,name="e%s"%self.filter)
        
        if "NUM" not in av.colnames:
            narr= np.nansum(np.invert(np.isnan(tab2array(tab, find_col_names(tab, "RA")))), axis=1)
            av.add_column(Column(narr, name="NUM"))

            if num_thresh>0:
                av.remove_rows( av["NUM"]<num_thresh)
        return av
    
    def build_meta(self,catalogues):
        """
        Not happy with this yet
        """
        meta=catalogues[0].meta
        base=Table(None, meta=meta)
        return base

        

class CascadeMatch(GenericMatch):
    """
    A simple advancement on "Generic Matching" where the number
    of columns are not preserved. At the end of each sub match, the 
    table is left justified, to reduce the total number of columns needed.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs, method="Cascade Matching")

    def match(self, catalogues, **kwargs):
        """
        Match a list of catalogues with RA and DEC columns

        Parameters
        ----------
        See `GenericMatch.match`

        Returns
        -------
        output : `astropy.table.Table`
            A left aligned catalogue of all the matched values
        """
        catalogues=self.init_catalogues(catalogues)
        if CAT_NUM in self.col_names: self.col_names.remove(CAT_NUM)
        base=self.build_meta(catalogues)

        for n,cat in enumerate(catalogues,1):
            self.load.msg="matching: %d"%n
            tmp=self._match(base,cat, join_type="or")
            tmp.rename_columns( tmp.colnames, ["%s_%d"%(name,n) for name in tmp.colnames] )
            base=h_cascade((base, tmp), col_names=self.col_names)
        base=fill_nan(base)
        return base

class DitherMatch(GenericMatch):
    """
    The same as Generic Matching

    """
    
    def __init__(self, catalogues, pfile=None, method="DitherMatch"):
        super(self,DitherMatch).__init__(catalogues, pfile)

    def match(self, **kwargs):
        """
        """
        return None

class BandMatch(GenericMatch):
    # filter flag for kwargs.
    FILTER = "fltr"
    THRESHOLD = "threshold"

    # match methods
    _FIRST = "first"
    _LAST = "last"
    _BOOT_STRAP = "bootstrap"

    # warning messages
    _WRONG_THRESHOLD = (
        "Threshold values must be scalar or list with length 1 less than the "
        "catalogue list. The final element is being ignored.\n")


    def __init__(self, **kwargs):
        if self.FILTER in kwargs:
            if not isinstance(kwargs[self.FILTER], list):
                warn("{} input should be a list, "
                     "there may be unexpected behaviour\n", self.FILTER)

        if self.THRESHOLD in kwargs:
            if isinstance(kwargs[self.THRESHOLD],list):
                kwargs[self.THRESHOLD]=np.array(kwargs[self.THRESHOLD])

        super().__init__(**kwargs, method="Band Matching")

    def order_catalogues(self,catalogues):
        """
        Reorder catalogue list into increasing wavelength size
        This only works for JWST bands. Unrecognised filters will be left unchanged.
        The function should also set the self.filter variable if possible

        Parameter
        ---------
        catalogues : list
            List of `astropy.table.Table` with meta keys FILTER

        Returns
        -------
        catalogues : list
            The same list reordered
        """

        status=-1
        _ii=None
        sorters = [
            ## META in JWST filters
            lambda t: list(starbug2.filters.keys()).index(t.meta.get(_FILTER)),
            ## colnames in JWST filters
            lambda t: list(starbug2.filters.keys()).index(
                (set(t.col_names) & set(starbug2.filters.keys())).pop()),
            ## META in self.filters
            lambda t: self.filter.index( t.meta.get(_FILTER)),
            ## colnames in JWST filters
            lambda t: self.filter.index(
                (set(t.col_names) & set(self.filter)).pop())
        ]

        for n, fn in enumerate(sorters):
            try:
                catalogues.sort(key=fn)
                _ii = map(fn,catalogues)
                status = n
                break
            except Exception as e:
                pass

        if status < 0:
            p_error(
                "Unable to reorder catalogues, leaving input order"
                " untouched.\n")
        ## JWST filters
        elif status <= 1 and (_ii is not None):
            self.filter=[list(starbug2.filters.keys())[i] for i in _ii]

        self.load=Loading(sum(len(c) for c in catalogues[1:]))

        return catalogues

    def jwst_order(self,catalogues):
        pass

    def match(self, catalogues, method="first", **kwargs):
        """
        Given a list of catalogues, it will reorder them into increasing
        wavelength or to match the fltr= keyword in the initializer.
        The matching then uses the shortest wavelength available position.
        I.e If F115W, F444W, F770W are input, the F115W centroid positions will
        be taken as "correct". If a source is not resolved in this band, the
        next most astrometric ally accurate position is taken, i.e. F444W

        :param catalogues: List of `astropy.table.Table` objects containing
        the meta item "FILTER=XXX"
        :type catalogues: List of `astropy.table.Table`
        :param method: Centroid method
            "first" -   Use the position corresponding to the earliest
                        appearance of the source
            "last"  -   Use the position corresponding to the latest
                        appearance of the source
            "bootsrap"- ..
            "average" - ..
        :type method: str
        :param kwargs:
        :return: Matched catalogue
        :rtype: astropy.Table
        """
        catalogues = self.order_catalogues(catalogues)

        if isinstance(self.filter,list) and len(self.filter)==len(catalogues):
            printf("Bands: %s\n"%', '.join(self.filter))
        else:
            printf("Bands: Unknown\n")

        if type(self.threshold.value) in (list,np.ndarray):
            if len(self.threshold) != (len(catalogues) - 1):
                warn(self._WRONG_THRESHOLD)
                self.threshold = self.threshold[:-1]
        else: 
            self.threshold = (
                np.full(len(catalogues) - 1, self.threshold)*u.arcsec)

        printf("Thresholds: %s\n"%", ".join(
            ["%g\""%g for g in self.threshold.value]))

        if self.col_names is None:
            self.col_names = [
                "RA", "DEC", "flag", "NUM",
                *self.filter, *["e%s" % f for f in self.filter]]

        printf("Columns: %s\n"%", ".join(self.col_names))

        if method not in (self._FIRST, self._LAST, self._BOOT_STRAP):
            method = self._FIRST

        #########
        # Begin #
        #########
        
        base=self.build_meta(catalogues)
        _threshold=self.threshold.copy()
        for n,tab in enumerate(catalogues):
            ## Temporarily recast threshold
            self.threshold = _threshold[n-1]
            self.load.msg="%s (%g\")" % (self.filter[n], self.threshold.value)
            col_names = [
                name for name in self.col_names if name in tab.col_names]

            tmp = self._match(base,tab, join_type="or")

            col_names.remove("RA")
            col_names.remove("DEC")
            base = fill_nan(hstack((base, tmp[col_names])))
            base.rename_columns(
                col_names, ["%s_%d"%(name,n+1) for name in col_names])

            if "RA" not in base.col_names:
                base = fill_nan(hstack((tmp["RA", "DEC"], base)))
            elif method == self._FIRST:
                _mask=np.logical_and( np.isnan(base["RA"]), tmp["RA"]!=np.nan)
                base["RA"][_mask] = tmp["RA"][_mask]
                base["DEC"][_mask] = tmp["DEC"][_mask]
            elif method == self._LAST:
                _mask = ~np.isnan(tmp["RA"])
                base["RA"][_mask] = tmp["RA"][_mask]
                base["DEC"][_mask] = tmp["DEC"][_mask]
            elif method == self._BOOT_STRAP:
                _mask = ~np.isnan(tmp["RA"])
                base.rename_columns(("RA","DEC"), ("_RA_%d"%n, "_DEC_%d"%n))
                base = hstack( (base, tmp[["RA","DEC"]]) )

        self.threshold=_threshold # Set threshold back at the end
        ####################
        # Fix column names #
        for name in self.col_names:
            all_cols=find_col_names(base, name)
            if len(all_cols)==1:
                base.rename_column(all_cols.pop(), name)

        ################################
        # Finalise NUM and flag column #
        tmp=self.finish_matching( base, colnames=["NUM", "flag"] )
        base.remove_columns(
            (*find_col_names(base, "NUM"), *find_col_names(base, "flag")))
        base.add_column(tmp["NUM"], index=2)
        base.add_column(tmp["flag"],index=3)

        return base


def band_match(catalogues, col_names=("RA","DEC")):
    """
    Given a list of catalogues (with filter names in the metadata), match
    them in order of decreasing astrometric accuracy. If F115W, F444W, F770W
    are input, the F115W centroid positions will be taken as "correct". If a
    source is not resolved in this band, the next most astrometric ally
    accurate position is taken, i.e. F444W

    :param catalogues:
    :param col_names:
    :return:
    """

    ### ORDER the tables into the correct order (increasing wavelength)
    tables = np.full( len(starbug2.filters), None)
    mask = np.full(len(starbug2.filters), False)
    for tab in catalogues:
        if _FILTER in tab.meta.keys():
            if tab.meta[_FILTER] in starbug2.filters:
                ii = list(starbug2.filters.keys()).index(tab.meta[_FILTER])
                tables[ii] = tab
                mask[ii] = True
            else:
                p_error(
                    "Unknown filter '%s' (skipping)..\n" % tab.meta[_FILTER])
        elif _tmp := set(starbug2.filters.keys()) & set(tab.col_names):
            ii = list(starbug2.filters.keys()).index(_tmp.pop())
            tables[ii] = tab
            mask[ii] = True
        else:
            p_error("Cannot find 'FILTER' in table meta (skipping)..\n")

    # document bands
    s = "Bands: "
    for filter_string, tab in zip(starbug2.filters.keys(),tables):
        if tab: s+="%5s "% filter_string
    puts(s)

    ### Match in increasing wavelength order
    base = Table(None)
    load = Loading(
        sum([len(t) for t in tables[mask][1:]]), "matching", res=100)
    for filter_string, tab in zip(starbug2.filters.keys(), tables):
        if not tab:
            continue

        ## removing empty magnitude rows
        tab.remove_rows(np.isnan(tab[filter_string]))
        load.msg = "matching:%s" % filter_string
        _col_names = list(name for name in tab.col_names if name in col_names)
        if not len(base): 
            tmp = tab[_col_names].copy()
        else:
            idx, d2d, _ = GenericMatch._match(base,tab)
            tmp = Table(
                np.full( (len(base),len(_col_names)), np.nan),
                names = _col_names)
            
            ###################################
            # Hard coding separations for now #
            separation = 0.06
            f_id = list(starbug2.filters.keys()).index(filter_string)
            if f_id >= list(starbug2.filters.keys()).index("F277W"):
                separation = 0.10
            if f_id >= list(starbug2.filters.keys()).index("F560W"):
                separation = 0.15
            if f_id >= list(starbug2.filters.keys()).index("F1000W"):
                separation = 0.20
            if f_id >= list(starbug2.filters.keys()).index("F1500W"):
                separation = 0.25

            for ii, (src, IDX, sep) in enumerate(zip(tab, idx, d2d)):
                load.msg = "matching:%s(%.2g\")" % (filter_string, separation)
                load()
                load.show()

                if ((sep <= separation * u.arcsec)
                        and (sep == min(d2d[idx == IDX]))):
                    for name in _col_names: tmp[IDX][name] = src[name]
                else:
                    tmp.add_row(src[_col_names])

        tmp.rename_column("flag", "flag_%s"%filter_string)
        base = hstack((
            base, tmp[[filter_string, "e%s"%filter_string,
            "flag_%s"%filter_string]]
        ))
        base = Table(base, dtype=[float]*len(base.col_names)).filled(np.nan)

        ### Only keep the most astronomically correct position
        if "RA" not in base.col_names:
            base = hstack((tmp[["RA", "DEC"]], base))
        else:
            _mask = np.logical_and( np.isnan(base["RA"]), tmp["RA"] != np.nan)
            base["RA"][_mask] = tmp["RA"][_mask]
            base["DEC"][_mask] = tmp["DEC"][_mask]

    ## Sort out flags
    flag = np.zeros(len(base),dtype=np.uint16)
    for f_col in find_col_names(base, "flag"):
        flag |= base[f_col].value.astype(np.uint16)
        base.remove_column(f_col)
    base.add_column(flag,name="flag")
    
    return base.filled(np.nan)


class ExactValueMatch(GenericMatch):
    """
    Match catalogues based on a single value in one of their columns
    The normal use case is matching sources based on their *Catalogue_Number*
    value.
    """

    def __init__(self, value=CAT_NUM, **kwargs):
        """
        setup method.

        :param value: Column name to take exact values from
        :type value: str
        :param kwargs:
        """
        self.value = value
        super().__init__(**kwargs, method="Exact Value Matching")

        if "colnames" in kwargs: 
            p_error("Colnames not implemented in %s\n" % self.method)

    def __str__(self):
        s=[ "%s:" % self.method,
            "Value: \"%s\"" % self.value,
            "Colnames: %s" % self.col_names,
            ]
            
        return "\n".join(s)

    def _match(self, base, cat):
        """
         The low level matching function.

        :param base: Table onto which to match *cat*
        :type base: astropy.Table
        :param cat: Table to match to *base*
        :type cat:  astropy.Table
        :return: A new catalogue, it is a reordered version of *cat*, in the
            correct sorting to be h-stacked with *base*
        :rtype:  astropy.Table
        """
        
        tmp = Table(
            np.full((len(base), len(cat.col_names)), np.nan),
            names=cat.col_names, dtype=cat.dtype, masked=True)
        for col in tmp.columns.values():
            col.mask |= True

        if not len(base):
            return vstack([tmp,cat])

        for src in cat:
            if self.verbose:
                self.load()
                self.load.show()
            ii = np.where(base[self.value] == src[self.value])[0]
            if len(ii):
                tmp[ii] = src
            else:
                tmp.add_row(src)
        return tmp

    def match(self, catalogues, **kwargs):
        """
        Core matching function

        :param catalogues: List of astropy Tables to match together
        :type catalogues: list of astropy.Table
        :param kwargs:
        :return: Full matched catalogue
        :rtype: astropy.Table
        """

        catalogues = self.init_catalogues(catalogues)
        base = self.build_meta(catalogues)

        if self.value not in self.col_names:
            p_error("Exact value '%s' not in column names.\n" % self.value)
            return None

        for n, cat in enumerate(catalogues, 1):
            self.load.msg = "matching: %d"%n
            tmp = self._match(base,cat)
            tmp.rename_columns(
                tmp.col_names, ["%s_%d" % (name, n) for name in tmp.col_names])
            base = hstack([base,tmp])

            if n > 1:
                ii = (base[self.value].mask
                      & ~base["%s_%d" % (self.value, n)].mask)
                base["%s" % self.value][ii] = (
                    base["%s_%d" % (self.value, n)][ii])
                base.remove_column("%s_%d" % (self.value, n))

            else: base.rename_column("%s_1" % self.value, self.value)
            
        return fill_nan(base)

def sort_exposures(catalogues):
    """
     Given a list of catalogue files, this will return the fitsHDULists as a
     series of nested dictionaries sorted by:
    >   BAND
    >   OBSERVATION ID
    >   VISIT ID
    >   DETECTOR                -- These two have been switched
    >   DITHER (EXPOSURE)       -- These two have been switched

    :param catalogues: the catalogues to sort exposures of.
    :return: a dictionary of sorted catalogues
    """
    out = {}
    for cat in catalogues:
        info = exp_info(cat)
        
        if info[_FILTER] not in out.keys():
            out[info[_FILTER]] = {}

        if info[_OBS] not in out[info[_FILTER]].keys():
            out[info[_FILTER]][info[_OBS]] = {}

        if info[_VISIT] not in out[info[_FILTER]][info[_OBS]].keys():
            out[info[_FILTER]][info[_OBS]][info[_VISIT]] = {}

        if (info[_DETECTOR] not in
                out[info[_FILTER]][info[_OBS]][info[_VISIT]].keys()):
            out[info[_FILTER]][
                info[_OBS]][info[_VISIT]][info[_DETECTOR]] = []
        out[info[_FILTER]][
            info[_OBS]][info[_VISIT]][info[_DETECTOR]].append(cat)
    return out


def parse_mask(string, table):
    """
    Parse an commandline mask string to be passed into a matching routine
    Example: --mask=F444W!=nan

    :param string: Raw mask sting to be parsed
    :type string: str
    :param table: Table to work on
    :type table: astropy.table.Table
    :return: Boolean mask array to index into a table or array
    :rtype: np.ndarray
    """
    mask=None
    
    for col_name in table.colnames: string=string.replace(
        col_name, "table[\"%s\"]" % col_name)

    try:
        mask = eval(string)
        if not isinstance(mask, np.ndarray):
            raise Exception
    except NameError as e:
        p_error("Unable to create mask: %s\n" % repr(e))
    except Exception as e:
        p_error(repr(e))

    return mask


def exp_info(hdu_list):
    """
    Get the exposure information about a hdu list
    INPUT:  HDUList or ImageHDU or BinTableHDU
    RETURN: dictionary of relevant information:
            >   EXPOSURE, DETECTOR, FILTER
    """
    info={  _FILTER : None,
            _OBS : 0,
            _VISIT : 0,
            _EXPOSURE : 0,
            _DETECTOR : None
            }

    if type(hdu_list) in (fits.ImageHDU, fits.BinTableHDU):
        hdu_list=fits.HDUList(hdu_list)

    for hdu in hdu_list:
        for key in info:
            if key in hdu.header:
                info[key] = hdu.header[key]
    return info