function JobSubmission(func,varargin)
%JobSubmission creates and submits dos batch or pbs scripts to relevant machine to execute certain functions.
%
% Usage: JobSubmission('RunMainScript','flag',value)
%    or: JobSubmission RunMainScript flag value
%
%e.g.
%Within Matlab
%             JobSubmission('RunMainScript ','--grid-para-buffer',1,'--dem-para-filename','/g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz','--grav-grad-para-avail',true);
%
%Compiled
%             JobSubmission RunMainScript --grid-para-buffer 1 --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz --grav-grad-para-avail true
%
%Available options:
%--grid-para-buffer <value>                             e.g. --grid-para-buffer 1
%--grid-para-buffer2 <value>                            e.g. --grid-para-buffer2 0.5
%--grid_para-step <value>                               e.g. --grid_para-step 0.5
%--grid-para-filtersize <value>                         e.g. --grid-para-filtersize 15
%--grid-para-filterradius <value>                       e.g. --grid-para-filterradius 10
%--grid-para-minlong <value>                            e.g. --grid-para-minlong 153
%--grid-para-maxlong <value>                            e.g. --grid-para-maxlong 154
%--grid-para-minlat <value>                             e.g. --grid-para-minlat -29
%--grid-para-maxlat <value>                             e.g. --grid-para-maxlat -28
%--dem-para-filename <path_to_file>                     e.g. --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz
%--dem-para-num-cols <value>                            e.g. --dem-para-num-cols 4861
%--dem-para-num-rows <value>                            e.g. --dem-para-num-rows 3181
%--grav-para-filename <path_to_file>                    e.g. --grav-para-filename /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW.mat
%--grav-para-filename1 <path_to_file>                   e.g. --grav-para-filename1 /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW_1.mat
%--grav-para-typeb <value>                              e.g. --grav-para-typeb 1
%--grav-para-grav-faye-typeb <value>                    e.g. --grav-para-grav-faye -typeb 3
%--grav-para-altimetry-weighting <value>                e.g. --grav-para-altimetry-weighting 1
%--grav-grad-para-filename <path_to_file>               e.g. --grav-grad-para-filename /g/data/dg9/nd2979/Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat
%--grav-grad-para-typeb <value>                         e.g. --grav-grad-para-typeb 0.00001
%--grav-grad-para-avail <logical>                       e.g. --grav-grad-para-avail true
%--cov-para-compute-empircal-cov-dec <value>            e.g. --cov-para-compute-empircal-cov-dec 3
%--cov-para-fit-empircal-cov <type>                     e.g. --cov-para-fit-empircal-cov auto
%--cov-para-fitempiricalcovnsearch <values>             e.g. --cov-para-fitempiricalcovnsearch 21600,1,21600
%--cov-para-fitempiricalcovmsearch <values>             e.g. --cov-para-fitempiricalcovmsearch 200,20,300
%--cov-para-n <value>                                   e.g. --cov-para-n 10800
%--cov-para-m <value>                                   e.g. --cov-para-m 200
%--cov-para-width <value>                               e.g. --cov-para-width 3
%--cov-para-res <value>                                 e.g. --cov-para-res 0.00833333333
%--cov-para-cov-computed_tilewise <logical>             e.g. --cov-para-cov-computed_tilewise true
%--cov-para-airbornedataonly <logical>                  e.g. --cov-para-airbornedataonly false
%--cov-para-covplot <logical>                           e.g. --cov-para-covplot false
%--topo-para-corr <logical>                             e.g. --topo-para-corr true
%--topo-para-topoplot <logical>                         e.g. --topo-para-topoplot false
%--topo-para-density <value>                            e.g. --topo-para-density 2.67
%--topo-para-depth <value>                              e.g. --topo-para-depth 0
%--topo-para-rad <value>                                e.g. --topo-para-rad 1
%--topo-para-rtm <values>                               e.g. --topo-para-rtm 50,10,300
%--ggm-para-filename <path_to_file>                     e.g. --ggm-para-filename /g/data/dg9/nd2979/Data/GGM/GOCE_For_Gridded_Int.mat
%--coast-para-filename <path_to_file>                   e.g. --coast-para-filename /g/data/dg9/nd2979/Data/COASTLINE/CoastAus.mat
%--levelling-para-lev-eval <logical>                    e.g. --levelling-para-lev-eval true
%--levelling-para-filename <path_to_file>               e.g. --levelling-para-filename /g/data/dg9/nd2979/Data/GPS_LEVELLING/Lev_NSW_NG.mat
%--levelling-para-plot-stats <logical>                  e.g. --levelling-para-plot-stats false
%--levelling-para-compare-to-existing-model <logical>   e.g. --levelling-para-compare-to-existing-model true
%--levelling-para-existing-model <path_to_file>         e.g. --levelling-para-existing-model /g/data/dg9/nd2979/Data/EXISTING_GEOID_MODELS/AGQG20221120.mat
%--levelling-para-max-diff <value>                      e.g. --levelling-para-max-diff 0.15
%--output-para-grids-name <path_to_folder>              e.g. --output-para-grids-name /g/data/dg9/nd2979/outputs/GridsNENSWgg2degTile/
%--output-para-tiles-dir-name <path_to_folder>          e.g. --output-para-tiles-dir-name /g/data/dh8/outputs/ResidualTilesNENSWgg2degTile/
%--output-para-plot-grids <logical>                     e.g. --output-para-plot-grids false
%--output-para-plotsfolder <path_to_folder>             e.g. --output-para-plotsfolder /g/data/dh8/outputs/plots/22-Nov-2024NENSWgg2degTile
%--keepawake <logical>                                  e.g. --keepawake true
%--memtype <string>                                     e.g. --memtype normal
%--nsubsets <value>                                     e.g. --nsubsets 8
%--walltime <string>                                    e.g. --walltime 48:00:00
%--mem <string>                                         e.g. --mem 512GB
%--jobfs <string>                                       e.g. --jobfs 400GB
%--ncpus <value>                                        e.g. --ncpus
%
%Geoscience Australia. Neda Darbeheshti on 24/11/2024
%Modified on 07/05/2025 by Justy Siwabessy
%

% Default parameters
execpath = '/g/data/dg9/gravityLibrary/executables/';
proj = 'dg9';
walltime = '24:00:00';
mem = '512GB';
jobfs = '400GB';
ncpus = '48';
memtype = 'hugemem';
ntasks = 1;
isnsubsets = false;
issubset_size = false;
dxdy = [];
str = [];
GRID_PARA.STEP=0.5; % The step size. This must be less than buffer2 to avoid gaps in the final grid.

%check for first input argument if --help
if strncmp(func,'--help',6)
    helptext
    return
end

%check for input argument
for i=1:nargin-1
    if ischar(varargin{i})
        if strncmp(varargin{i},'--grid-para-buffer',18)
            GRID_PARA.buffer = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid-para-buffer2',19)
            GRID_PARA.buffer2 = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid_para-step',16)
            GRID_PARA.STEP = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid-para-filtersize',22)
            GRID_PARA.filterSize = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid-para-filterradius',24)
            GRID_PARA.filterRadius = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid-para-minlong',19)
            GRID_PARA.MINLONG = str2num(varargin{i+1});
            %str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid-para-maxlong',19)
            GRID_PARA.MAXLONG = str2num(varargin{i+1});
            %str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid-para-minlat',18)
            GRID_PARA.MINLAT = str2num(varargin{i+1});
            %str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grid-para-maxlat',18)
            GRID_PARA.MAXLAT = str2num(varargin{i+1});
            %str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--dem-para-filename',19)
            DEM_PARA.filename = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--dem-para-num-cols',19)
            DEM_PARA.num_cols = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--dem-para-num-rows',19)
            DEM_PARA.num_rows = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-para-filename',20)
            GRAV_PARA.filename = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-para-filename1',21)
            GRAV_PARA.filename1 = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-para-typeb',17)
            GRAV_PARA.TypeB = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-para-grav-faye-typeb',27)
            GRAV_PARA.Grav_Faye_TypeB = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-para-altimetry-weighting',21)
            GRAV_PARA.altimetry_weighting = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-grad-para-filename',25)
            GRAV_GRAD_PARA.filename = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-grad-para-typeb',22)
            GRAV_GRAD_PARA.TypeB = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--grav-grad-para-avail',22)
            GRAV_GRAD_PARA.avail = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-compute-empircal-cov-dec',35)
            COV_PARA.Compute_Empircal_COV_Dec = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-fit-empircal-cov',27)
            COV_PARA.Fit_Empircal_COV = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-fitempiricalcovnsearch',33)
            COV_PARA.FitEmpiricalCOVNSearch = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-fitempiricalcovmsearch',33)
            COV_PARA.FitEmpiricalCOVMSearch = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-n',12)
            COV_PARA.N = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-m',12)
            COV_PARA.M = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-width',16)
            COV_PARA.width = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-res',14)
            COV_PARA.res = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-cov-computed_tilewise',32)
            COV_PARA.COV_COMPUTED_Tilewise = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-airbornedataonly',27)
            COV_PARA.Airbornedataonly = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--cov-para-covplot',18)
            COV_PARA.COVPlot = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--topo-para-corr',16)
            Topo_PARA.Corr = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--topo-para-topoplot',20)
            Topo_PARA.TopoPlot = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--topo-para-density',19)
            Topo_PARA.Density = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--topo-para-depth',17)
            Topo_PARA.Depth = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--topo-para-rad',15)
            Topo_PARA.Rad = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--topo-para-rtm',15)
            Topo_PARA.RTM = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--ggm-para-filename',19)
            GGM_PARA.filename = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--coast-para-filename',21)
            COAST_PARA.filename = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--levelling-para-lev-eval',25)
            LEVELLING_PARA.Lev_eval = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--levelling-para-filename',25)
            LEVELLING_PARA.filename = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--levelling-para-plot-stats',27)
            LEVELLING_PARA.Plot_Stats = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--levelling-para-compare-to-existing-model',42)
            LEVELLING_PARA.Compare_To_Existing_Model = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--levelling-para-existing-model',31)
            LEVELLING_PARA.Existing_Model = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--levelling-para-max-diff',25)
            LEVELLING_PARA.max_diff = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--output-para-grids-name',24)
            OUTPUT_PARA.Grids_name = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--output-para-tiles-dir-name',28)
            OUTPUT_PARA.Tiles_dir_name = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--output-para-plot-grids',24)
            OUTPUT_PARA.PLOT_GRIDS = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--output-para-plotsfolder',25)
            OUTPUT_PARA.plotsFolder = varargin{i+1};
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strncmp(varargin{i},'--output-para-polygonLon',24)
            OUTPUT_PARA.polygonLon = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--output-para-polygonLat',24)
            OUTPUT_PARA.polygonLat = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strncmp(varargin{i},'--keepawake',11)
            keepawake = str2num(varargin{i+1});
            str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--executables-path',18)
            execpath = varargin{i+1};
            % str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--proj',6)
            proj = varargin{i+1};
            % str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--memtype',9)
            memtype = varargin{i+1};
            % str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--walltime',10)
            walltime = varargin{i+1};
            % str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--mem',5)
            mem = varargin{i+1};
            % str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--jobfs',7)
            jobfs = varargin{i+1};
            % str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--ncpus',7)
            ncpus = varargin{i+1};
            % str = [str ' ' varargin{i} ' ' varargin{i+1} ''];
        elseif strncmp(varargin{i},'--nsubsets',10)
            ntasks = str2num(varargin{i+1});
            if ~ntasks
                disp('Must be > 0');
                return
            end
            isnsubsets = true;
            issubset_size = false;
        elseif strncmp(varargin{i},'--subset-size',13)
            dxdy = str2num(varargin{i+1});
            if ~dxdy
                disp('Must be > 0');
                return
            elseif length(dxdy) > 2
                disp('Length must be < 3')
                return
            end
            issubset_size = true;
            isnsubsets = false;
      elseif strncmp(varargin{i},'--help',6)
            helptext
            return
        end
    end
end

% terminate if both --nsubsets and --subset-size have been selected.
if isnsubsets & issubset_size
    disp('Cannot be both --nsubsets and --subset-size')
    return
end

if isnsubsets
    subset_boundaries = divide_area_flexible(GRID_PARA.MINLONG, GRID_PARA.MAXLONG, GRID_PARA.MINLAT, GRID_PARA.MAXLAT, ntasks);
elseif issubset_size
    subset_boundaries = divide_area_by_cell_size(GRID_PARA.MINLONG, GRID_PARA.MAXLONG, GRID_PARA.MINLAT, GRID_PARA.MAXLAT, dxdy(1), dxdy(end));
    ntasks = length(subset_boundaries);
    ncpus = num2str(dxdy(1) / GRID_PARA.STEP * dxdy(end) / GRID_PARA.STEP);
end
c4pbs = ['walltime=' walltime ',' 'mem=' mem ',' 'jobfs=' jobfs ',' 'ncpus=' ncpus]; %pbs setting
joblist = {};

for i=1:1:ntasks
    nowinsec = round(rem(now,1)*24*60*60*1000);
    boundarystr = [' --grid-para-minlong ' num2str(subset_boundaries{i}(1))  ' --grid-para-maxlong ' num2str(subset_boundaries{i}(2))  ' --grid-para-minlat ' num2str(subset_boundaries{i}(3))  ' --grid-para-maxlat ' num2str(subset_boundaries{i}(4))];
    if strncmpi(computer,'pcwin',5)
        jobfname = [func '_' num2str(nowinsec) '.bat'];
        fidjob = fopen(jobfname,'w');
        % fprintf(fidjob,'%s\n','set PATH=%PATH%;C:\W10DEV\app\neda');
        fprintf(fidjob,'%s\n',['set PATH=%PATH%;' execpath]);
        fprintf(fidjob,'%s\n','prompt $N:\$G');
        fprintf(fidjob,'%s\n','');
        fprintf(fidjob,'%s\n',[func '' str '' boundarystr]);
        fclose(fidjob);
        % dos([jobfname ' &']);
    else
        jobfname = [func '_' num2str(nowinsec) '.job'];
        fidjob = fopen(jobfname,'w');
        fprintf(fidjob,'%s\n','#!/bin/bash');
        fprintf(fidjob,'%s\n',['#PBS -P ' proj]);
        fprintf(fidjob,'%s\n',['#PBS -q ' memtype]);
        fprintf(fidjob,'%s\n',['#PBS -l ' c4pbs]);
        fprintf(fidjob,'%s\n','#PBS -l wd');
        fprintf(fidjob,'%s\n','#PBS -l storage=gdata/dg9');
        fprintf(fidjob,'%s\n','');
        % fprintf(fidjob,'%s\n','export PATH="/home/547/jxs547/executables:$PATH"'); %always check for the one used and loaded in raijin
        fprintf(fidjob,'%s\n',['export PATH="' execpath ':$PATH"']); %always check for the one used and loaded in raijin
        fprintf(fidjob,'%s\n','module load matlab'); %always check for the one used and loaded in raijin
        % fprintf(fidjob,'%s\n',['~/executables/' func '' str '' boundarystr]);
        fprintf(fidjob,'%s\n',[func '' str '' boundarystr]);
        fclose(fidjob);
        % unix(['qsub ' jobfname]);
    end
    joblist(i,:) = cellstr(jobfname);
end

if exist('joblist','var');
    for i=1:1:size(joblist,1)
        if strncmpi(computer,'pcwin',5)
            % dos(['type ' joblist{i}]);
            dos([joblist{i} ' &']);
        else
            % unix(['cat ' joblist{i}]);
            unix(['qsub ' joblist{i}]);
        end
    end
    return
end

function subset_boundaries = divide_area_flexible(min_lon, max_lon, min_lat, max_lat, num_divisions)
%DIVIDE_AREA_FLEXIBLE Divides a rectangular area into a flexible number of subsets.
%
%   subset_boundaries = DIVIDE_AREA_FLEXIBLE(MIN_LON, MAX_LON, MIN_LAT, MAX_LAT, NUM_DIVISIONS)
%
%   Inputs:
%       MIN_LON        Minimum longitude of the area.
%       MAX_LON        Maximum longitude of the area.
%       MIN_LAT        Minimum latitude of the area.
%       MAX_LAT        Maximum latitude of the area.
%       NUM_DIVISIONS  The desired total number of subsets to divide the
%                      area into. The function will try to divide as evenly
%                      as possible into a grid-like structure.
%
%   Output:
%       SUBSET_BOUNDARIES A cell array where each cell contains a 1x4
%                         vector defining the boundaries of a subset in the
%                         order [min_lon, max_lon, min_lat, max_lat].

subset_boundaries = {};

if num_divisions <= 0
    error('Number of divisions must be a positive integer.');
end

% Try to find the closest integer factors for rows and columns
factors = factor(num_divisions);
n_rows = 1;
n_cols = 1;

if ~isempty(factors)
    % Start with factors closest to the square root
    sqrt_n = sqrt(num_divisions);
    best_diff = Inf;
    best_r = 1;
    best_c = num_divisions;

    for r = 1:num_divisions
        if mod(num_divisions, r) == 0
            c = num_divisions / r;
            diff = abs(r - c);
            if diff < best_diff
                best_diff = diff;
                best_r = r;
                best_c = c;
            elseif diff == best_diff && r < best_r % Prefer more rows if difference is the same
                best_r = r;
                best_c = c;
            end
        end
    end
    n_rows = best_r;
    n_cols = best_c;
else
    % For prime numbers, divide along one dimension
    n_rows = 1;
    n_cols = num_divisions;
end

delta_lon = (max_lon - min_lon) / n_cols;
delta_lat = (max_lat - min_lat) / n_rows;

count = 0;
for i = 0:n_rows-1
    for j = 0:n_cols-1
        count = count + 1;
        current_min_lon = min_lon + j * delta_lon;
        current_max_lon = min_lon + (j + 1) * delta_lon;
        current_min_lat = min_lat + i * delta_lat;
        current_max_lat = min_lat + (i + 1) * delta_lat;
        subset_boundaries{count} = [current_min_lon, current_max_lon, current_min_lat, current_max_lat];
    end
end

function subset_boundaries = divide_area_by_cell_size(min_lon, max_lon, min_lat, max_lat, cell_width_lon, cell_height_lat)
%DIVIDE_AREA_BY_CELL_SIZE Divides a rectangular area into subsets based on desired cell dimensions.
%
%   subset_boundaries = DIVIDE_AREA_BY_CELL_SIZE(MIN_LON, MAX_LON, MIN_LAT, MAX_LAT, CELL_WIDTH_LON, CELL_HEIGHT_LAT)
%
%   Inputs:
%       MIN_LON           Minimum longitude of the overall area (in degrees).
%       MAX_LON           Maximum longitude of the overall area (in degrees).
%       MIN_LAT           Minimum latitude of the overall area (in degrees).
%       MAX_LAT           Maximum latitude of the overall area (in degrees).
%       CELL_WIDTH_LON    Desired longitude span (width) for each individual subset/cell (in degrees).
%                         Must be a positive value.
%       CELL_HEIGHT_LAT   Desired latitude span (height) for each individual subset/cell (in degrees).
%                         Must be a positive value.
%
%   Output:
%       SUBSET_BOUNDARIES A cell array where each cell contains a 1x4
%                         vector defining the boundaries of a subset in the
%                         order [min_lon, max_lon, min_lat, max_lat].
%                         The subsets are ordered row by row (bottom to top),
%                         then column by column (left to right).

% Validate input cell dimensions
if cell_width_lon <= 0 || cell_height_lat <= 0
    error('Cell width (CELL_WIDTH_LON) and height (CELL_HEIGHT_LAT) must be positive values.');
end

% Calculate total span of the area
total_lon_span = max_lon - min_lon;
total_lat_span = max_lat - min_lat;

% Calculate the total number of columns and rows needed to cover the entire area.
% Use 'ceil' to ensure that even if the span is not perfectly divisible by the
% cell size, enough cells are created to cover the whole area.
n_cols_total = ceil(total_lon_span / cell_width_lon);
n_rows_total = ceil(total_lat_span / cell_height_lat);

% Pre-allocate the cell array for efficiency
subset_boundaries = cell(1, n_cols_total * n_rows_total);
count = 0;

% Iterate through rows (latitude divisions) from bottom to top
for i = 0:n_rows_total-1
    % Iterate through columns (longitude divisions) from left to right
    for j = 0:n_cols_total-1
        count = count + 1;

        % Calculate the min and max longitude for the current cell
        current_min_lon = min_lon + j * cell_width_lon;
        current_max_lon = min_lon + (j + 1) * cell_width_lon;

        % Calculate the min and max latitude for the current cell
        current_min_lat = min_lat + i * cell_height_lat;
        current_max_lat = min_lat + (i + 1) * cell_height_lat;

        % Important: Cap the maximum longitude/latitude of the last cells
        % to ensure they don't extend beyond the original overall boundary.
        % This prevents going slightly over max_lon/max_lat due to ceil.
        current_max_lon = min(current_max_lon, max_lon);
        current_max_lat = min(current_max_lat, max_lat);

        % Store the boundaries
        subset_boundaries{count} = [current_min_lon, current_max_lon, current_min_lat, current_max_lat];
    end % End of inner (j) for-loop
end % End of outer (i) for-loop

function helptext
str={'JobSubmission creates and submits dos batch or pbs scripts to relevant machine to execute certain functions.'
' '
' Usage: JobSubmission(''RunMainScript'',''flag'',value)'
'    or: JobSubmission RunMainScript flag value'
' '
'e.g.'
'Within Matlab'
'             JobSubmission(''RunMainScript '',''--grid-para-buffer'',1,''--dem-para-filename'',''/g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz'',''--grav-grad-para-avail'',true);'
' '
'Compiled'
'            JobSubmission --grid-para-buffer 1 --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz --grav-grad-para-avail true'
' '
'Available options:'
'--grid-para-buffer <value>                             e.g. --grid-para-buffer 1'
'--grid-para-buffer2 <value>                            e.g. --grid-para-buffer2 0.5'
'--grid_para-step <value>                               e.g. --grid_para-step 0.5'
'--grid-para-filtersize <value>                         e.g. --grid-para-filtersize 15'
'--grid-para-filterradius <value>                       e.g. --grid-para-filterradius 10'
'--grid-para-minlong <value>                            e.g. --grid-para-minlong 153'
'--grid-para-maxlong <value>                            e.g. --grid-para-maxlong 154'
'--grid-para-minlat <value>                             e.g. --grid-para-minlat -29'
'--grid-para-maxlat <value>                             e.g. --grid-para-maxlat -28'
'--dem-para-filename <path_to_file>                     e.g. --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz'
'--dem-para-num-cols <value>                            e.g. --dem-para-num-cols 4861'
'--dem-para-num-rows <value>                            e.g. --dem-para-num-rows 3181'
'--grav-para-filename <path_to_file>                    e.g. --grav-para-filename /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW.mat'
'--grav-para-filename1 <path_to_file>                   e.g. --grav-para-filename1 /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW_1.mat'
'--grav-para-typeb <value>                              e.g. --grav-para-typeb 1'
'--grav-para-grav-faye-typeb <value>                    e.g. --grav-para-grav-faye -typeb 3'
'--grav-para-altimetry-weighting <value>                e.g. --grav-para-altimetry-weighting 1'
'--grav-grad-para-filename <path_to_file>               e.g. --grav-grad-para-filename /g/data/dg9/nd2979/Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat'
'--grav-grad-para-typeb <value>                         e.g. --grav-grad-para-typeb 0.00001'
'--grav-grad-para-avail <logical>                       e.g. --grav-grad-para-avail true'
'--cov-para-compute-empircal-cov-dec <value>            e.g. --cov-para-compute-empircal-cov-dec 3'
'--cov-para-fit-empircal-cov <type>                     e.g. --cov-para-fit-empircal-cov auto'
'--cov-para-fitempiricalcovnsearch <values>             e.g. --cov-para-fitempiricalcovnsearch 21600,1,21600'
'--cov-para-fitempiricalcovmsearch <values>             e.g. --cov-para-fitempiricalcovmsearch 200,20,300'
'--cov-para-n <value>                                   e.g. --cov-para-n 10800'
'--cov-para-m <value>                                   e.g. --cov-para-m 200'
'--cov-para-width <value>                               e.g. --cov-para-width 3'
'--cov-para-res <value>                                 e.g. --cov-para-res 0.00833333333'
'--cov-para-cov-computed_tilewise <logical>             e.g. --cov-para-cov-computed_tilewise true'
'--cov-para-airbornedataonly <logical>                  e.g. --cov-para-airbornedataonly false'
'--cov-para-covplot <logical>                           e.g. --cov-para-covplot false'
'--topo-para-corr <logical>                             e.g. --topo-para-corr true'
'--topo-para-topoplot <logical>                         e.g. --topo-para-topoplot false'
'--topo-para-density <value>                            e.g. --topo-para-density 2.67'
'--topo-para-depth <value>                              e.g. --topo-para-depth 0'
'--topo-para-rad <value>                                e.g. --topo-para-rad 1'
'--topo-para-rtm <values>                               e.g. --topo-para-rtm 50,10,300'
'--ggm-para-filename <path_to_file>                     e.g. --ggm-para-filename /g/data/dg9/nd2979/Data/GGM/GOCE_For_Gridded_Int.mat'
'--coast-para-filename <path_to_file>                   e.g. --coast-para-filename /g/data/dg9/nd2979/Data/COASTLINE/CoastAus.mat'
'--levelling-para-lev-eval <logical>                    e.g. --levelling-para-lev-eval true'
'--levelling-para-filename <path_to_file>               e.g. --levelling-para-filename /g/data/dg9/nd2979/Data/GPS_LEVELLING/Lev_NSW_NG.mat'
'--levelling-para-plot-stats <logical>                  e.g. --levelling-para-plot-stats false'
'--levelling-para-compare-to-existing-model <logical>   e.g. --levelling-para-compare-to-existing-model true'
'--levelling-para-existing-model <path_to_file>         e.g. --levelling-para-existing-model /g/data/dg9/nd2979/Data/EXISTING_GEOID_MODELS/AGQG20221120.mat'
'--levelling-para-max-diff <value>                      e.g. --levelling-para-max-diff 0.15'
'--output-para-grids-name <path_to_folder>              e.g. --output-para-grids-name /g/data/dg9/nd2979/outputs/GridsNENSWgg2degTile/'
'--output-para-tiles-dir-name <path_to_folder>          e.g. --output-para-tiles-dir-name /g/data/dh8/outputs/ResidualTilesNENSWgg2degTile/'
'--output-para-plot-grids <logical>                     e.g. --output-para-plot-grids false'
'--output-para-plotsfolder <path_to_folder>             e.g. --output-para-plotsfolder /g/data/dh8/outputs/plots/22-Nov-2024NENSWgg2degTile'
'--keepawake <logical>                                  e.g. --keepawake true'
'--memtype <string>                                     e.g. --memtype normal'
'--nsubsets <value>                                     e.g. --nsubsets 8'
'--walltime <string>                                    e.g. --walltime 48:00:00'
'--mem <string>                                         e.g. --mem 512GB'
'--jobfs <string>                                       e.g. --jobfs 400GB'
'--ncpus <value>                                        e.g. --ncpus'
' '
'Geoscience Australia. Neda Darbeheshti on 24/11/2024'
'Modified on 07/05/2025 by Justy Siwabessy'};

for i=1:length(str)
    disp(str{i})
end
return
