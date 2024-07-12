TREESKIP = 100
prefix = "penguinsrelate/pop_penguins_popsize_chr{CHR}"
anc = f"{prefix}.anc"
mut = f"{prefix}.mut"
dist = f"{prefix}.dist"
coal = "penguinsrelate/pop_penguins_popsize.coal"
# CHRS = [1, 2, 3, 4, 5]
CHRS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
# CHRS = [18]
ms = ['1.5e-8']
nsampless = [100]
# tCutoffs = [None, int(1e6), int(1e5), int(1e4)] 
# tCutoffs = [int(1e4)]
tCutoffs = [None]


# ------------ decide which trees to sample ----------------

bps = anc.replace('anc','bps').replace('penguinsrelate','data') 

rule get_bps:
  input:
    expand(bps, CHR=CHRS, m=ms)

checkpoint get_bp:
  input:
    anc, mut
  output:
    bps 
  run:
    print('getting tree indices')
    ixs_start=[]
    ixs_end=[]
    with open(input[0], "r") as f:
      for i, line in enumerate(f): #each line is a tree, see https://myersgroup.github.io/relate/getting_started.html#Output
        if i==1: 
          n = int(line.split()[1]) #number of trees on this chromosome
          trees = [i for i in range(0,n+1,TREESKIP)] #which trees to sample
        if i > 1 and i-2 in trees:
          ixs_start.append(int(line.split(':')[0])) #index of first snp in sampled tree
        if i > 2 and i-3 in trees: 
          ixs_end.append(int(line.split(':')[0])-1) #index of last snp in sampled tree
    print('choose',len(ixs_start),'trees')
    print('getting start and stop basepairs')
    bps_start = []
    bps_end = []
    with open(input[1],"r") as f:
      for i,line in enumerate(f):
        if i>0 and int(line.split(';')[0]) in ixs_start:
          bps_start.append(int(line.split(';')[1])) #position of first snp in sampled tree
        if i>0 and int(line.split(';')[0]) in ixs_end:
          bps_end.append(int(line.split(';')[1])) #position of last snp in sampled tree
    print('writing to file')
    with open(output[0], 'w') as out:
      for start,end in zip(bps_start,bps_end):
        out.write(str(start) + ' ' + str(end) + '\n')
        
# ------------ sample branch lengths at a particular location -------------

tree = bps.replace('.bps','_{start}-{stop}bps_{nsamples}nsamples.newick')

def input_func(name):

  def input_files(wildcards):
    filenames = []
    for CHR in CHRS:
        for m in ms:
          infile = checkpoints.get_bp.get(CHR=CHR, m=m).output[0]
          with open(infile,'r') as f:
            for line in f:
              i,j = line.strip().split(' ')
              d = {'{CHR}': CHR, '{m}': m, '{start}': i, '{stop}': j}
              string = name
              for i,j in d.items():
                string = string.replace(i,str(j))
              filenames.append(string)
    return expand(filenames, nsamples=nsampless)

  return input_files

rule sample_trees:
  input:
    input_func(tree) 

ruleorder: sample_tree > get_bp

rule sample_tree:
    input:
        expand(f"{prefix}.{{end}}", end=["anc", "mut", "dist"], allow_missing=True),
        coal
    output:
        tree 
    params:
        prefix_in = f"{prefix}",
        prefix_out = tree.replace('.newick','')
    threads: 1
    group: 'sample_tree'
    resources: 
        runtime = 15
    shell:
        '''
        module load gcc/13.2.0
        echo "DEBUG: input={input}, params={params}"  
        ~/scratch/relate/scripts/SampleBranchLengths/SampleBranchLengths.sh \
                     -i {params.prefix_in} \
                     --dist {input[2]} \
                     --coal {input[3]} \
                     -o {params.prefix_out}  \
                     -m 1.5e-8 \
                     --format n \
                     --num_samples {wildcards.nsamples} \
                     --first_bp {wildcards.start} \
                     --last_bp {wildcards.stop} \
                     --seed 1 
        '''
# snakemake sample_trees --groups sample_tree=sample_tree --group-components sample_tree=80 --jobs 1

# ------------ shared time matrices and coalescence times -------------

shared_times = tree.replace('.newick','_sts.npy')
coal_times = tree.replace('.newick','_cts.npy')

rule times:
  input:
    input_func(shared_times),
    input_func(coal_times)

rule time:
  input:
    tree
  output:
    shared_times,
    coal_times
  threads: 1
  resources: 
        runtime = 15
  group: 'time'
  run:
    from tsconvert import from_newick
    from utils import get_shared_times
    import numpy as np
    from tqdm import tqdm

    stss = []
    ctss = []
    with open(input[0], mode='r') as f:
      next(f) #skip header
      for line in tqdm(f, total=int(wildcards.nsamples)): #for each tree sampled at a locus

        # import
        string = line.split()[4] #extract newick string only (Relate adds some info beforehand)
        ts = from_newick(string) #convert to tskit "tree sequence" (only one tree)
        tree = ts.first() #the only tree

        # shared times
        samples = [int(ts.node(node).metadata['name']) for node in ts.samples()] #get index of each sample in list we gave to relate
        sample_order = np.argsort(samples) #get indices to put in ascending order
        ordered_samples = [ts.samples()[i] for i in sample_order] #order samples as in relate
        sts = get_shared_times(tree, ordered_samples) #get shared times between all pairs of samples, with rows and columns ordered as in relate
        stss.append(sts)

        #coalescence times
        cts = sorted([tree.time(i) for i in tree.nodes() if not tree.is_sample(i)]) #coalescence times, in ascending order
        ctss.append(cts)

    np.save(output[0], np.array(stss))
    np.save(output[1], np.array(ctss))
    
# ------------- get locations ---------------------

rule locations:
    input:
        population="penguinsrelate/output_data.txt",
        site_data="penguinsrelate/site_data.csv"
    output:
        "penguinsrelate/matched_data.npy"
    run:
        import csv
        import numpy as np
        from pyproj import CRS, Transformer

        def convert_projected_to_lonlat(x, y, source_epsg, dest_epsg=4326):
            source_proj = pyproj.Proj(f'epsg:{source_epsg}')
            dest_proj = pyproj.Proj(f'epsg:{dest_epsg}')
            lon, lat = pyproj.transform(source_proj, dest_proj, x, y)
            return lat, lon

        population_data = []
        with open(input.population, 'r') as pop_file:
            for line in pop_file:
                sample, population = line.strip().split()
                population_data.append((sample, population))

        site_data = {}
        source_epsg = 3031  
        with open(input.site_data, 'r') as site_file:
            reader = csv.DictReader(site_file)
            for row in reader:
                long, lat = float(row['LONG']), float(row['LAT'])
                lon, lat = convert_projected_to_lonlat(long, lat, source_epsg)
                site_data[row['SITE_CODE']] = (lon, lat)

        matched_data = []
        for sample, population in population_data:
            if population in site_data:
                lon, lat = site_data[population]
                matched_data.append([lon, lat])

        match_array = np.array(matched_data)
        matched_array = np.repeat(match_array, 2, axis=0)
        np.save(output[0], matched_array)

        csv_file_path = "penguinsrelate/matched_data.csv"
        with open(csv_file_path, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            for row in matched_array:
                writer.writerow(row)


# ------------- process shared times ----------

processed_shared_times = shared_times.replace('_sts.npy','_sts_{tCutoff}tCutoff_{end}.npy')
ends = ['mc-invs','logdets','samples']

rule process_shared_time:
  input:
    shared_times
  output:
    expand(processed_shared_times, end=ends, allow_missing=True)
  threads: 1
  group: 'process_shared_time'
  resources: 
        runtime = 15
  run:
    # taming numpy
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    # load times at all trees
    import numpy as np
    stss = np.load(input[0])
    _,n,_ = stss.shape
    # process trees 
    from utils import chop_shared_times, center_shared_times
    stss_inv = []
    stss_logdet = []
    smplsss = []
    tCutoff = wildcards.tCutoff
    if tCutoff=='None': 
      tCutoff=None 
    else:
      tCutoff=float(tCutoff)
    for sts in stss:
      # chop
      sts_chopped, smpls = chop_shared_times(sts, tCutoff=tCutoff) #shared times and samples of each subtree
      sts_inv = []
      sts_logdet = []
      smplss = []
      # process subtrees
      for st,sm in zip(sts_chopped, smpls):
        stc = center_shared_times(st) #mean center
        stc_inv = np.linalg.pinv(stc) #invert
        stc_logdet = np.linalg.slogdet(st)[1] #log determinant
        sts_inv.append(stc_inv)
        sts_logdet.append(stc_logdet) 
        smplss.append(sm) #samples
      stss_inv.append(sts_inv)
      stss_logdet.append(sts_logdet)
      smplsss.append(smplss)
    # save
    np.save(output[0], np.array(stss_inv, dtype=object))
    np.save(output[1], np.array(stss_logdet, dtype=object))
    np.save(output[2], np.array(smplsss, dtype=object))
    
# ------------- process coalescence times ----------

processed_coal_times = coal_times.replace('_cts.npy','_cts_{tCutoff}tCutoff_{end}.npy')
ends = ['bts','lpcs']

rule process_coal_time:
  input:
    coal_times,
    coal
  output:
    expand(processed_coal_times, end=ends, allow_missing=True)
  threads: 1 
  group: 'process_coal_time'
  resources: 
        runtime = 15
  run:
    import os
    # taming numpy
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    import numpy as np
    from spacetrees import _log_coal_density
    # get variable Nee
    epochs = np.genfromtxt(input[1], skip_header=1, skip_footer=1) #time at which each epoch starts (and the final one ends)
    Nes = 0.5/np.genfromtxt(input[1], skip_header=2)[2:] #effective population size during each epoch
    # process coal times
    ctss = np.load(input[0]) #coalescence times in ascending order, for each tree
    btss = []
    lpcs = []
    tCutoff = wildcards.tCutoff
    if tCutoff=='None': 
      tCutoff=None 
    else:
      tCutoff=float(tCutoff)
    for cts in ctss: 
      # get branching times in ascending order
      T = cts[-1] #TMRCA (to be replaced with tCutoff)
      if tCutoff is not None:
        if tCutoff < T:
          T = tCutoff
      bts = T - np.flip(cts) #branching times, in ascending order
      bts = bts[bts>0] #remove branching times at or before T
      bts = np.append(bts,T) #append total time as last item
      btss.append(bts)
      # get probability of coalescence times under panmictic coalescent with variable Ne
      lpc = _log_coal_density(times=cts, Nes=Nes, epochs=epochs, tCutoff=tCutoff)
      lpcs.append(lpc)
    # save
    np.save(output[0], np.array(btss, dtype=object))
    np.save(output[1], np.array(lpcs, dtype=object))

# ------------- composite dispersal rates -------------

locations = "penguinsrelate/matched_data.npy"

composite_dispersal_rate = processed_shared_times.replace('chr{CHR}_','').replace('_{start}-{stop}bps','').replace('_sts_{tCutoff}tCutoff_{end}','_{tCutoff}tCutoff_mle-dispersal')

rule composite_dispersal_rates:
  input:
    expand(composite_dispersal_rate, nsamples=nsampless, m=ms, tCutoff=tCutoffs)

def input_func_dispersal(name, ends):
  def input_files(wildcards):
    filenames = []
    for CHR in CHRS:
        bpfile = checkpoints.get_bp.get(CHR=CHR, **wildcards).output[0]
        with open(bpfile, 'r') as f:
          for line in f:
            start, stop = line.strip().split(' ')
            d = {'{CHR}': CHR, '{start}': start, '{stop}': stop}
            string = name
            for i, j in d.items():
              string = string.replace(i, str(j))
            filenames.append(string)
    return expand(filenames, end=ends, allow_missing=True)
  return input_files

rule composite_dispersal_rate:
  input:
    stss_mc_inv = input_func_dispersal(processed_shared_times, ['mc-invs']),
    stss_logdet = input_func_dispersal(processed_shared_times, ['logdets']),
    smplss = input_func_dispersal(processed_shared_times, ['samples']),
    btss = input_func_dispersal(processed_coal_times, ['bts']),
    lpcss = input_func_dispersal(processed_coal_times, ['lpcs']),
    locs = expand(locations, CHR=CHRS, allow_missing=True)[0]
  output:
    composite_dispersal_rate
  threads: 80
  group: 'composite_dispersal_rate'
  resources:
    runtime = 300
  run:
    import os
    # taming numpy
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    import numpy as np
    from spacetrees import mle_dispersal, _sds_rho_to_sigma
    from tqdm import tqdm
    # load locations
    locations = np.load(input.locs)
    print("Size of locations array:", locations.shape[0])
    # subsample for testing
    L = len(input.stss_mc_inv)
    M = int(wildcards.nsamples)
    # L = 10 #number of loci
    # M = 10 #number of trees per locus
    # mean centered and inverted shared time matrices
    print('\nloading inverted shared times matrices')
    stss_mc_inv = []
    for f in tqdm(input.stss_mc_inv[:L]):
      sts_mc_inv = np.load(f, allow_pickle=True)[:M]
      stss_mc_inv.append(sts_mc_inv)
    # log determinants of mean centered shared time matrices    
    print('\nloading log determinants of shared times matrices')
    stss_logdet = []
    for f in tqdm(input.stss_logdet[:L]):
      sts_logdet = np.load(f, allow_pickle=True)[:M]
      stss_logdet.append(sts_logdet) 
    # subtree samples    
    print('\nloading samples of shared times matrices')
    smplss = []
    for f in tqdm(input.smplss[:L]):
      smpls = np.load(f, allow_pickle=True)[:M]
      smplss.append(smpls) 
    # branching times
    print('\nloading branching times')
    btss = []
    for f in tqdm(input.btss[:L]):
      bts = np.load(f, allow_pickle=True)[:M] 
      btss.append(bts)
    # log probability of coalescent times    
    print('\nloading log probability of coalescence times')
    lpcss = []
    for f in tqdm(input.lpcss[:L]):
      lpcs = np.load(f, allow_pickle=True)[:M] 
      lpcss.append(lpcs)
    # function for updates
    def callbackF(x):
    #  print('{0: 3.6f}   {1: 3.6f}   {2: 3.6f}'.format(x[0], x[1], x[2])) #if important=False
      print('{0: 3.6f}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}'.format(x[0], x[1], x[2], x[3]))
    # find parameter estimates
    print('\nestimating dispersal rate')
    mle = mle_dispersal(locations=locations, shared_times_inverted=stss_mc_inv, log_det_shared_times=stss_logdet, samples=smplss, 
                        sigma0=None, phi0=None, # make educated guess based on first tree at each locus
                        # sigma0=_sds_rho_to_sigma(0.07, 0.06, 0.5), phi0=5e-5, # guess from tCutoff=1e6 mle
                        callbackF=callbackF, 
                        important=True, branching_times=btss, logpcoals=lpcss)
    print('\n', mle)
    np.save(output[0], mle)
    
# ------------- ancestor locations using chopped trees -------------

ancestor_locations = processed_shared_times.replace('_sts','').replace('{end}','anc-locs')

def input_func_locs(name):

  def input_files(wildcards):
    filenames = []
    for CHR in CHRS:
        for m in ms:
          infile = checkpoints.get_bp.get(CHR=CHR, m=m).output[0]
          with open(infile,'r') as f:
            for line in f:
              i,j = line.strip().split(' ')
              d = {'{CHR}': CHR, '{m}': m, '{start}': i, '{stop}': j}
              string = name
              for i,j in d.items():
                string = string.replace(i,str(j))
              filenames.append(string)
    return expand(filenames, nsamples=nsampless, tCutoff=tCutoffs)

  return input_files

rule locate_ancestors:
  input:
    input_func_locs(ancestor_locations)

ruleorder: process_shared_time > locate_ancestor 
ruleorder: process_coal_time > locate_ancestor 

rule locate_ancestor:
  input:
    stss = shared_times,
    locs = locations,
    mle = composite_dispersal_rate,
    btss = processed_coal_times.replace('{end}','bts'),
    lpcs = processed_coal_times.replace('{end}','lpcs'),
  output:
    ancestor_locations 
  threads: 1 #get seg fault if >1
  group: 'locate_ancestor'
  resources: 
        runtime = 200
  run:
    import os
    # taming numpy
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    
    import numpy as np
    from spacetrees import _sds_rho_to_sigma, _log_birth_density, locate_ancestors
    from utils import chop_shared_times

    # sample locations
    locations = np.load(input.locs)
    n = len(locations)

    # who and when
    ancestor_samples = range(n) #which samples to find ancestors of
    tCutoff = wildcards.tCutoff
    if tCutoff == 'None':
      tCutoff = None
      max_time = 1e6
    else:
      tCutoff = float(tCutoff)
      max_time = tCutoff  
    ancestor_times = np.logspace(1,np.log10(max_time),10) #times to find ancestors

    # chop trees
    stss_chopped = []
    samples = []
    for sts in np.load(input.stss):
      sts_chopped, smpls = chop_shared_times(sts, tCutoff=tCutoff) #shared times and sample indices in each subtree
      stss_chopped.append(sts_chopped)
      samples.append(smpls)
    
    #dispersal rate
    mle = np.load(input.mle) #mle dispersal rate and branching rate
    sigma = _sds_rho_to_sigma([mle[0], mle[1], mle[2]]) #as covariance matrix

    # importance weights
    btss = np.load(input.btss, allow_pickle=True) #birth times
    phi = mle[-1] #mle branching rate
    lbds = np.array([_log_birth_density(bts, phi, n) for bts in btss]) #log probability densities of birth times
    lpcs = np.load(input.lpcs, allow_pickle=True) #log probability densities of coalescence times
    log_weights = lbds - lpcs #log importance weights
 
    # subsample for testing
    M = int(wildcards.nsamples)
    #M = 100 #number of trees

    # locate 
    ancestor_locations = locate_ancestors(ancestor_samples, ancestor_times, 
                                          stss_chopped[:M], samples[:M], locations, log_weights[:M], sigma)
    np.save(output[0], ancestor_locations)
# taking about 16m each with 1 thread
# snakemake locate_ancestors --profile slurm --groups locate_ancestor=locate --group-components locate=80 --jobs 11

# ------------- ancestor locations using full trees -------------

ancestor_locations_full = processed_shared_times.replace('_sts','').replace('{end}','anc-locs_full-trees')

rule locate_ancestors_full:
  input:
    input_func_locs(ancestor_locations_full)

ruleorder: process_shared_time > locate_ancestor_full
ruleorder: process_coal_time > locate_ancestor_full

rule locate_ancestor_full:
  input:
    stss = shared_times,
    locs = locations,
    mle = composite_dispersal_rate,
    btss = processed_coal_times.replace('{end}','bts'),
    lpcs = processed_coal_times.replace('{end}','lpcs'),
  output:
    ancestor_locations_full 
  threads: 1 #get seg fault if >1
  run:
    import os
    # taming numpy
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    
    import numpy as np
    from spacetrees import _sds_rho_to_sigma, _log_birth_density, locate_ancestors
    from utils import chop_shared_times

    # sample locations
    locations = np.load(input.locs)
    n = len(locations)

    # who and when
    ancestor_samples = range(n) #which samples to find ancestors of
    tCutoff = wildcards.tCutoff #determines which dispersal rate we use and how far back we locate ancestors
    if tCutoff == 'None':
      max_time = 1e6
    else:
      max_time = float(tCutoff)  
    #ancestor_times = np.logspace(1,np.log10(max_time),10) #times to find ancestors
    ancestor_times = np.linspace(1e3,max_time,10)

    # chop trees
    stss_chopped = []
    samples = []
    for sts in np.load(input.stss):
      sts_chopped, smpls = chop_shared_times(sts, tCutoff=None) #shared times and sample indices in each subtree (here just 1 subtree per tree)
      stss_chopped.append(sts_chopped)
      samples.append(smpls)
    
    #dispersal rate
    mle = np.load(input.mle) #mle dispersal rate and branching rate
    sigma = _sds_rho_to_sigma([mle[0], mle[1], mle[2]]) #as covariance matrix

    # importance weights
    btss = np.load(input.btss, allow_pickle=True) #birth times
    phi = mle[-1] #mle branching rate
    lbds = np.array([_log_birth_density(bts, phi, n) for bts in btss]) #log probability densities of birth times
    lpcs = np.load(input.lpcs, allow_pickle=True) #log probability densities of coalescence times
    log_weights = lbds - lpcs #log importance weights
 
    # locate 
    ancestor_locations = locate_ancestors(ancestor_samples, ancestor_times, 
                                          stss_chopped, samples, locations, log_weights, sigma)
    np.save(output[0], ancestor_locations)
# taking about 16m each with 1 thread
# snakemake locate_ancestors_full --profile slurm --groups locate_ancestor_full=locate --group-components locate=80 --jobs 11

