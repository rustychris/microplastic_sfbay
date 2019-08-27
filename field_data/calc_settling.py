import numpy as np
import logging as log

# extra handling:
#  Fiber Bundle: dimensions are probably not great, since it's going to be very porous 

density={}
# density["Rubber"]=1.522 # https://www.aqua-calc.com/page/density-table/substance/rubber-coma-and-blank-manufactured
# density["Rubber"]=1.2 # https://www.engineeringtoolbox.com/density-solids-d_1265.html
# These are primarily fragments
density["Rubber"]=1.15 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5664766/, referencing a Fed. Hwy. Admin. report

density["Polyethylene"]=0.94 # middle of range https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Cotton"]=1.50 # confirm? https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Cellulose acetate"]=1.30 # sheet, https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Polypropylene"]=0.9 # https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Cellulosic"]=1.50 # https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Polyvinyl chloride"]=1.39 # # https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Acrylic"]=1.19
density["Organic natural material"]=density["Cellulosic"] # using same
density["Glass"]=2.60 # https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Polyethylene co-vinyl acetate"]=0.94 # https://www.sigmaaldrich.com/catalog/substance/polyethylenecovinylacetate123452493778811?lang=en&region=US&attrlist=Density
# PVB is typically a film in safety glass.
density["Polyvinyl butyral"]=1.08 # http://www2.dupont.com/Building_Innovations/zh_CN/assets/downloads/SGPintro_E.pdf
density["Nylon"]=1.15 # https://www.engineeringtoolbox.com/density-solids-d_1265.html
density["Cellulosic (natural-based)"]=density['Cellulosic']
density["Anthropogenic (cellulosic)"]=density['Cellulosic'] # fibers
density["Polystyrene"]=1.03 # unless expanded... https://www.engineeringtoolbox.com/density-solids-d_1265.html 
density["Polyethylene terephthalate"]=1.38 # https://en.wikipedia.org/wiki/Polyethylene_terephthalate
density["Polyacrylamide"]=1.30 # https://polymerdatabase.com/polymers/polyacrylamide.html
density["Polyester"]=1.39 # https://nptel.ac.in/courses/116102026/36
density["Polyamide"]=1.14 # https://www.azom.com/article.aspx?ArticleID=477
density["Wool"]=1.29 # https://www.researchgate.net/publication/47929535_Measurement_of_density_and_medullation_in_wool
density["Paint"]=1.4 # is this alkyd binder? https://pubs.acs.org/doi/pdf/10.1021/es501757s
# paints were mostly fragments.

density["Methyl vinyl ether copolymers"]=1.05 # http://polymerdatabase.com/polymers/polymethylvinylether.html

density["Polyurethane"]=1.1 # as a foam, or solid?  some of both. https://en.wikipedia.org/wiki/Polyurethane
# density depends on details of production, but 1.1 is a reasonable central value.

density["Acrylonitrile butadiene styrene"]=1.07 # https://en.wikipedia.org/wiki/Acrylonitrile_butadiene_styrene#cite_note-MATBASE-1

# 0.94 vs 0.90 for polyethylene and polypropylene, so split the difference
density["Polyethylene/polypropylene copolymer"]=0.92

# density["Anthropogenic (unknown base)"]=X # primarily fibers
# density["Anthropogenic (synthetic)"]=0.5 # primarily foams.
# density["Inorganic natural material"]=X

density["Polyethylene co-acrylic acid"]=0.96 # https://www.sigmaaldrich.com/catalog/product/aldrich/181048?lang=en&region=US
density["Polytetrafluoroethylene"]=2.2 # https://en.wikipedia.org/wiki/Polytetrafluoroethylene

density["Polyvinyl acetate"]=1.19 # https://en.wikipedia.org/wiki/Polyvinyl_acetate

density['Phenolic resin']=1.3 # https://www.sciencedirect.com/topics/chemical-engineering/phenolic-resins


# dev:
# mean_fiber_width=combined[ combined['Category_Final']=='Fiber']['Width.mm'].mean()
# module:
mean_fiber_width=0.09 # mm

def record_to_ws(rec):
    """
    returns settling velocity, i.e. positive=sinking
    """
    if rec['Category_Final']=='Foam':
        # assume these are buoyant foams
        rho_bulk=0.1
    elif rec['PlasticType_Final'] in density:
        rho_bulk=density[rec['PlasticType_Final']]
    else:
        # may be able to punt better later.
        return np.nan

    # estimate three dimensions, a,b,c in mm,
    # and d_equivalent in m.
    a=rec['Length.mm']
    b=rec['Width.mm']
    if a<b:
        # happens...
        a,b=b,a
        
    cat=rec['Category_Final']
    if cat=='Fiber':
        d_equi=b/1000. # as per Rise and Fall, mm->m
        c=b
        if np.isnan(a):
            log.info("Fabricating fiber length of 5mm")
            a=5.0
    else:
        if cat=='Film':
            c=0.05 # total guess
        elif cat=='Fiber Bundle':
            # may have a better way -- this makes it look like a film
            c=mean_fiber_width
        elif cat in ['Fragment','Sphere','Foam',np.nan]:
            c=b # maybe there is a better approximation?
        else:
            raise Exception("Unhandled category %s"%cat)
        d_equi=(a*b*c)**(1./3) / 1000. # mm -> m

    if np.isnan(d_equi) or d_equi==0.0:
        # could maybe get some info from the size class / sieve 
        return np.nan
        
    # Rise and Fall, equation 14
    g=9.8 # gravity
    rho_water=1.025 # average SF Bay, g/ml
    CSF=c/np.sqrt(a*b)
    if cat=='Sphere':
        P=6 # perfectly round
    else:
        # no information, so go for middle of the scale
        P=3
        
    nu=1.05e-6 # m2/s, 18degC.

    w_x=0.1 # m/s

    for i in range(1000):
        Re=w_x*d_equi/nu
        if rho_bulk>rho_water:
            if cat == 'Fiber': # what about Fiber Bundle??
                C_D= (4.7/np.sqrt(Re) + np.sqrt(CSF))
            else:
                C_D= 3./( CSF * Re**(1./3) )
        else:
            if cat == 'Fiber': # what about Fiber Bundle??
                C_D = (10/np.sqrt(Re) + np.sqrt(CSF))
            else:
                # how did they come up with this? gnar-gnar
                C_D = (20/Re + 10/np.sqrt(Re) + np.sqrt(1.195-CSF)) * (6/P)**(1-CSF)
        new_w_x=np.sqrt(4./3 * d_equi/C_D * np.abs(rho_bulk-rho_water)/rho_water * g)
        if np.abs(new_w_x - w_x)<1e-6:
            break
        w_x=new_w_x
    else:
        raise Exception("failed to converge")
        w_x=np.nan

    if rho_bulk < rho_water:
        w_x*=-1
    return w_x

# note that some particles have high surface tension and are essentially surface-bound
# even if the density is > 1.03
# d_equi for films is probably not good
# Weird that drag coefficient is different for rising vs falling particles.
