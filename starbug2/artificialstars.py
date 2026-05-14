import numpy as np
from photutils.datasets import make_model_image, make_random_models_table
from photutils.psf import FittableImageModel
from astropy.table import Table,hstack
from astropy.io import fits
from scipy.optimize import curve_fit

try: import matplotlib.pyplot as plt
except:
    import matplotlib; matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

from starbug2.utils import (
    printf, p_error, crop_hdu, get_mj_ysr2jy_scale_factor, warn)
from starbug2.matching import GenericMatch

class Artificial_StarsIII():
    """
    ast
    """
    def __init__(self, starbug, index=-1):
        ## Initials the starbug instance
        self.starbug=starbug
        _=self.starbug.image
        _=self.starbug.load_psf()

        self.psf=FittableImageModel(self.starbug.psf)
        self.index=index

    def __call__(self,*args,**kwargs): return self.auto_run(*args,**kwargs)

    def auto_run(self, ntests, stars_per_test=1, subimage_size=-1, mag_range=(18,27),
            loading_buffer=None, autosave=-1, skip_phot=0, skip_background=0):
        """
        The main entry point into the artificial star test
        This handles everything except the results compilation at the end,

        Parameters
        ----------
        ntests : int
            Number of tests to run

        stars_per_test : int
            Number of stars to inject per test

        subimage_size : int
            in prep.

        mag_range : tuple,list
            Length two list or tuple containing the magnitude range of
            injected stars. These will be uniformly sampled from within 
            this range.

        loading_buffer : numpy.ndarray
            Length 3 array of shared memory to increment a loading bar
            between multiple subprocesses..

        autosave : int
            Auto quick saving output frequency

        skip_phot : int
            If true then ignore the PSF phot step and use aperture fluxes instead

        skip_background : int
            If true then ignore the background substraction step
            
        Returns
        -------
        test_result : astropy.Table
            Full raw test results. Injected initial properties with measured values
        """

        test_result=Table(np.full((ntests*stars_per_test,8),np.nan), names=["x_0","y_0","mag","flux","x_det","y_det","flux_det", "status"])
        scalefactor= get_mj_ysr2jy_scale_factor(self.starbug.image)
        base_image=self.starbug._image.copy()
        base_shape=np.copy(self.starbug.image.shape)
        stars_per_test=int(stars_per_test)
        passed=0

        ZP = self.starbug.options.get("ZP_MAG") if self.starbug.options.get("ZP_MAG") else 0
        buffer=0

        if mag_range[0]-mag_range[1] >=0:
            warn("Detected magnitude range in wrong order, put bright limit first\n")
            return None

        if any(base_shape < subimage_size):
            subimage_size=min(base_shape)
            p_error("subimage size greater than image size, setting to 'safe' value %d.\n" % subimage_size)
            
        for test in range(1,int(ntests)+1):
            centre=0
            centre= (base_shape[0]*np.random.random(), base_shape[1]*np.random.random())

            #image=cropHDU( base_image.__deepcopy__(), (0,-1), (0,-1) )
            image=base_image.__deepcopy__()
            #image=self.create_subimage( base_image.__deepcopy__(), subimage_size, position=centre, hdu=self.st

            shape=image[self.starbug._n_hdu].shape

            sourcelist= make_random_models_table( stars_per_test, { "x_0":[buffer,shape[0]-buffer],
                                                                    "y_0":[buffer,shape[1]-buffer],
                                                                    "mag":mag_range}) 
            sourcelist.add_column( 10.0 ** ( (ZP-sourcelist["mag"])/2.5 ) , name="flux")
            sourcelist.remove_column("id")

            #image[self.starbug._nHDU].data*=0
            star_overlay=make_model_image( shape, self.psf, sourcelist,model_shape=self.psf.shape)/scalefactor
            image[self.starbug._n_hdu].data+=star_overlay
            self.starbug._image=image
            
            n=len(sourcelist)
            result=self.single_test(image, sourcelist, skip_phot=skip_phot, skip_background=skip_background)
            passed+=sum(result["status"])
            test_result[(test-1)*stars_per_test: test*stars_per_test]=result

            if loading_buffer is not None:
                loading_buffer[0]+=1
                loading_buffer[2]=int(100*passed/(test*stars_per_test))

            if autosave>0 and not test%autosave:
                test_result.write("sbast-autosave%d.tmp"%self.index, overwrite=True, format="fits")
            del image # is this neccessary?
        return test_result

    def single_test(self, image, contains, skip_phot=0, skip_background=0):
        """
        Conduct a single test on an image with a set of initial source properties

        Parameters
        ----------
        image : numpy.ndarray
            2D image array to conduct test on

        contains : table
            Table of initial source properties to be injected into the image.
            This table must contain the columns ("x_0","y_0","flux")

        skip_phot : int
            Skip the PSF phot routine

        skip_background : int
            Skip the background estimation and subtraction step

        Returns
        -------
        result : Table
            Table hoizontally stacked with the initial inputs and the detection and 
            photometric results. Plus column named "status", an integer flag as to 
            whether the source was detected or not.
        """
        NULL=0
        DETECT=1
        test_result=Table(np.full((len(contains),4),np.nan), names=["x_det","y_det","flux_det","status"])

        threshold=2
        if not self.starbug.detect(): #Run detection on the image
            det=self.starbug.detections
            for i, src in enumerate(contains): #Check for detection in output
                separations=np.sqrt( (src["x_0"]-det["xcentroid"])**2 + (src["y_0"]-det["ycentroid"])**2)
                best_match=np.argmin(separations)
                if separations[best_match]<threshold:
                    test_result["x_det"][i]=det["xcentroid"][best_match]
                    test_result["y_det"][i]=det["ycentroid"][best_match]
                    test_result["flux_det"][i]=det["flux"][best_match]
                    test_result["status"][i]=DETECT
                else: test_result["status"][i]=NULL

            if sum(test_result["status"]) and (skip_background or not self.starbug.bgd_estimate()):  # Run background estim if there were detections
                self.starbug.detections = test_result

                if not skip_phot and not self.starbug.photometry(): # Run PSF photometry on detected sources
                    self.starbug.psfcatalogue.rename_columns(("x_init","y_init","xydev"),("_x_init","_y_init","_xydev"))
                    matched=GenericMatch(threshold=threshold)([contains, self.starbug.psfcatalogue], cartesian=True)
                    test_result["flux_det"] = matched[:len(test_result)]["flux_2"]

        return hstack((contains,test_result))


    def create_subimage(self, image, size, position=(0,0), hdu=1, buffer=0):
        """
        probably to be deprecated
        """
        subimage=None
        imshape=00
        if size<=0: return image,0,0
        x_edge=0
        y_edge=0
        if any(imshape < size):
            size=min(imshape)
            p_error("subimage size greater than image size, setting to 'safe' value %d.\n" % size)

        x_edge = int(max( position[0]-(size/2), buffer ))
        y_edge = int(max( position[1]-(size/2), buffer ))
        x_end =  int(min( position[0]+(size/2), imshape[0]-buffer))
        y_end =  int(min( position[1]+(size/2), imshape[1]-buffer))

        return crop_hdu(image, x_limit=(x_edge, x_end), y_limit=(y_edge, y_end)), x_edge, y_edge

def get_completeness(test_result):
    """
    Compile the results into magnitude binned values of recovery fraction
    and flux error

    Parameters
    ----------
    test_result : table
        The output from auto_run

    Returns
    -------
    result : astropy Table
        Table containing percent completeness as a function of magnitude
    """

    bins = np.arange( np.floor(min(test_result["mag"])), np.ceil(max(test_result["mag"])), 0.1)
    percs= np.zeros(len(bins))
    errors=np.zeros(len(bins))
    offsets=np.zeros(len(bins))
    means =np.zeros(len(bins))
    
    ibins = np.digitize( test_result["mag"], bins=bins)
    for i in range(max(ibins)):
        binned=test_result[ (ibins==i) ]
        if binned: percs[i]=float(sum(binned["status"]))/len(binned)

        mag_inj= -2.5*np.log10( binned["flux"])
        mag_det= -2.5*np.log10( binned["flux_det"])
        errors[i]=np.nanstd( mag_inj-mag_det )
        means[i]=np.nanmean( mag_inj-mag_det )
        offsets[i]=np.nanmedian(binned["flux"]/binned["flux_det"])


    out=Table( [bins,percs,errors,offsets], names=("mag","rec","err","off"), dtype=(float,float,float,float))
    return out

def get_spatialcompleteness(test_result,image,res=10):
    """
    Produce an image array showing the spatially dependant recovery fraction 

    Parameters
    ----------
    test_result : table
        The output from auto_run

    image : numpy.ndarry
        2D image array to take the shape from 

    res : int
        The resolution of the spatial bins

    Returns
    -------
    percs : numpy.ndarray
        A 2D array the same shape as the image input, pixel values
        show the fraction of injected sources recovered in this bin
    """
    if not image: return None
    xbins=np.arange(min(test_result["x_0"]),max(test_result["x_0"]), int(res))
    ybins=np.arange(min(test_result["y_0"]),max(test_result["y_0"]), int(res))
    percs=np.zeros(image.shape)
    for xi in xbins[:-1]:
        xo=xi+res
        for yi in ybins[:-1]:
            yo=yi+res
            mask=(test_result["x_0"]>=xi) & (test_result["x_0"]<xo) &(test_result["y_0"]>=yi) & (test_result["y_0"]<yo) 
            binned=test_result[mask] 
            if len(binned): percs[int(xi):int(xo),int(yi):int(yo)]=float(sum(binned["status"])/len(binned))
    return percs

def estim_completeness_mag(ast):
    """
    Estimate the completenss level of the artificial star test
    
    Parameters
    ----------
    ast : astropy Table
        Output of Artificial_Stars.get_completeness, table must contain columns (mag, rec)

    Returns
    -------
    fit : list
        The fitting parameters to the logistic curve f(x)=l/(1+exp(-k(x-xo)))
        fit=[l,xo,k]

    complete : list
        Magnitude of 70% and 50% completeness 
    """
    fit=[None,None,None]
    compl=[None,None,None]
    fn_i=lambda y,l,k,xo: xo-(np.log((l/y)-1)/k)

    if len(set(ast.col_names) & set(("mag", "rec")))==2:
        try:
            fit,_=curve_fit(scurve, ast["mag"], ast["rec"], [1, -1,np.median(ast["mag"])])
            compl=(fn_i(0.9,*fit),fn_i(0.7,*fit),fn_i(0.5,*fit))
        except:
            warn("Unable to fit completeness fractions\n")
    else: p_error("Input table must have columns 'mag' and 'rec'\n")
    return fit,compl

def scurve(x,l,k,xo):
    """
    S-curve function to fit completeness results to

    f(x)=l/(1+exp(-k(x-xo)))

    Parameters
    ----------
    x : list
        Magnitude range to input into function

    l,xo,k : float
        Function parameters

    Returns
    -------
    f(x) : float
    """
    return l/(1+np.exp(-k*(x-xo)))

def compile_results(raw, image=None, plotast=None, fltr="m"):
    """
    Compile all the raw data into usable results

    Parameters
    ----------

    Returns
    -------
    """
    completeness=get_completeness(raw)
    _cfit,_compl=estim_completeness_mag(completeness)
    spatial_completeness = get_spatialcompleteness(raw, image, res=10)

    head={  "COMPLETE_FN":"F(x)=l/(1+exp(-k(x-xo)))",
            "l":_cfit[0], "k":_cfit[1], "xo":_cfit[2] }
    for i,frac in enumerate((90,70,50)):
        if _compl[i] and not np.isnan(_compl[i]):
            printf("-> complete to %d%%: %s=%.2f\n"%(frac,fltr,_compl[i]))
            head["COMPLETE %d%%"%frac]=_compl[i]

    results= fits.HDUList( [fits.PrimaryHDU( header=fits.Header(head)), 
                            fits.BinTableHDU(data=completeness, name="AST"),
                            fits.BinTableHDU(data=raw, name="RAW"),
                            fits.ImageHDU(data=spatial_completeness, name="CMP")])

    if plotast:
        fig,ax=plt.subplots(1,figsize=(3.5,3),dpi=300)
        ax.scatter(completeness["mag"],completeness["rec"], c='k', lw=0, s=8)
        ax.plot(completeness["mag"],scurve(completeness["mag"],*_cfit),c='g',label=r"$f(x)=\frac{%.2f}{1+e^{%.2f(x-%.2f)}}$"%(_cfit[0],-_cfit[1],_cfit[2]))
        ax.axvline(_compl[0], c="seagreen",ls='--', label=("90%%:%.2f"%_compl[0]),lw=0.75)
        ax.axvline(_compl[1], c="seagreen",ls='-.', label=("70%%:%.2f"%_compl[1]),lw=0.75)
        ax.axvline(_compl[2], c="seagreen",ls=':', label=("50%%:%.2f"%_compl[2]),lw=0.75)
        ax.scatter(_compl,(0.9,0.7,0.5),marker='*', c='teal', s=10)
        ax.tick_params(direction="in",top=True,right=True)
        ax.set_title("Artificial Star Test")
        ax.set_xlabel(fltr)
        ax.set_ylabel("Fraction Recovered")
        ax.set_yticks([0,.25,.5,.75,1])
        ax.legend(loc="lower left",frameon=False, fontsize=8)
        plt.tight_layout()
        fig.savefig(plotast,dpi=300)
        printf("--> %s\n"%plotast)

    return results
