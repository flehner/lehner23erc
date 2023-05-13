close all
clear all

% ---------------------------
% lehner23erc_fig3_fig5.m
% ---------------------------
% Code to produce Figs. 3 and 5 in Lehner and Deser (2023)
%
% Citation:
% Lehner, F., C. Deser (2023):
% Origin, importance, and predictive limits of internal climate variability
% Environmental Research: Climate, DOI: https://doi.org/10.1088/2752-5295/accf30
%
% Notes on code:
% - Only does calculations and plotting, not pre-processing
% - Requires pre-processed monthly mean data, specifically:
%   - Connected CMIP6 simulations from historical, ssp126, ssp245, ssp370, ssp585, bilinearly regridded to 2.5°x2.5°
%   - CESM1-LE data, bilinearly regridded to 2.5°x2.5°
%   - CESM1-DPLE data, bilinearly regridded to 2.5°x2.5°
% ---------------------------

% -- select which figures to plot:
fig3 = 1; % 0=no, 1=yes
fig5 = 1; % 0=no, 1=yes

load_dple = 1; % CESM1 DPLE is a big single file (13GB), so avoid if not needed


pathin      = '~/Dropbox/work/';
vars        = {'tas','pr','psl'};
comp        = 'Amon';
seasons     = {'DJF','JJA','annual','JJAS'};

% -- Scenarios
scen_cmip6 = {'ssp126','ssp245','ssp370','ssp585'};

% ------------------------------------------------------------------------------
% -- model selection
models_cmip6 = {'ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CAMS-CSM1-0',...
                'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1',...
                'CanESM5','CanESM5-CanOE','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L',...
                'FGOALS-g3','FIO-ESM-2-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G',...
                'HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G',...
                'MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
                'MRI-ESM2-0','NESM3','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
models_cmip6(ismember(models_cmip6,'GFDL-CM4')) = []; % missing ssp126 and ssp370
models_cmip6(ismember(models_cmip6,'NorESM2-LM')) = []; % pr: ssp585 r1i1p1f1 exists, but the historical r1i1p1f1 only starts in 1950 (?)
models_cmip6(ismember(models_cmip6,'FIO-ESM-2-0')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'HadGEM3-GC31-LL')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'NESM3')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'KACE-1-0-G')) = []; % piControl only 150 years

models_cmip6 = models_cmip6(1:3)

% ------------------------------------------------------------------------------
ensmem_cmip6  = ones(length(models_cmip6),1);
% ------------------------------------------------------------------------------

start_cmip5  = 1870;
ende_cmip5   = 2100;
start_cmip6  = 1850;
ende_cmip6   = 2100;
% -- general paramaters
start       = 1950;
ende        = 2099;
time        = start:ende;
refstart    = 2001;
refende     = 2020;
reflength   = refende-refstart+1;
lati        = [-88.75:2.5:88.75];
loni        = [1.25:2.5:358.75];
wl          = 10; % 1 10 % runnig mean length in years
% -- area --
filein  = [pathin 'cmip5-ng/area_g025.nc'];
area    = ncread([filein],'AREA');
filein  = [pathin 'cmip5-ng/landfrac_g025.nc'];
land    = ncread([filein],'LANDFRAC');
id = isnan(land);
land(id) = 0;
landarea = land.*area;


% -- location to plot --
regions     = {'global','uk','sahel','europe','southern_ocean',...
               'se_asia','us_southwest','nino34','north_america','kolkata',...
               'dallas','alaska','seattle','upper_colorado','sydney',...
               'NH','SH','arctic','AA','zurich',...
               'southern_europe','northern_europe','denver','anchorage','india',...
               'sahara','NoCal','SoCal','Mediterranean','Central_Europe',...
               'Northern_Europe','NCA_nw','UCRB','ireland','NCA_CONUS',...
               'Florida'};

% -- select region
for r =  [1] % 30 %[1 2 3 4 5 8 13 21 22] %1:length(regions) % <------------------------------------------- EDIT!
  region = regions{r}
  clear('lon2','var0','seas0')
  if strcmp(region,'NCA')==1
    lat         = [20 70];
    lon         = [-125 -65]+360;
  end
  if strcmp(region,'NCA_CONUS')==1
    lat         = [26.1 49.1];
    lon         = [-124.9+360 -67.1+360];
  end
  if strcmp(region,'upper_colorado')==1
    lat         = [38.75 41.25];
    lon         = [248.75 251.25 253.75];
      seas0     = 3;
      var0      = 2;
  end
  if strcmp(region,'Florida')==1
    lat         = [25 31];
    lon         = [-87 -80]+360;
  end
  if strcmp(region,'UCRB')==1
    lat         = [37 42];
    lon         = [-111 -105]+360;
  end
  if strcmp(region,'us_southwest')==1
    lat         = [32.5 42.1];
    lon         = [-125 -103]+360;
  end
  if strcmp(region,'NoCal')==1
    lat         = [37.1 42.1];
    lon         = [-125 -113]+360;
  end
  if strcmp(region,'SoCal')==1
    lat         = [32.5 37.1];
    lon         = [-125 -113]+360;
  end
  if strcmp(region,'europe')==1
    lat         = [36.4 70.1];
    lon         = [0 23.4];
    lon2        = [350 360];
  end
  if strcmp(region,'southern_europe')==1
    lat         = [36.4 50.1];
    lon         = [0 23.4];
    lon2        = [350 360];
      % seas0     = 2;
      % var0      = 1;
  end
  if strcmp(region,'northern_europe')==1
    lat         = [50.1 70.1];
    lon         = [0 23.4];
    lon2        = [350 360];
  end
  if strcmp(region,'se_asia')==1
    lat         = [7.9 30.1];
    lon         = [93.2 122.1];
  end
  if strcmp(region,'india')==1
    lat         = [6.2 30.1];
    lon         = [69.9 87.3];
      seas0     = 2;
      var0      = 2;
  end
  if strcmp(region,'nino34')==1
    lat         = [-5 5];
    lon         = [-170 -120]+360;
      seas0        = 1;
      var0        = 1;
  end
  if strcmp(region,'uk')==1
    lat         = [50.2:2.5:58.8];
    lon         = [-10.9:2.5:0]+360;
      seas0       = 3;
      var0        = 1;
  end
  if strcmp(region,'ireland')==1
    lat         = [52.2:2.5:55.8];
    lon         = [-10.9:2.5:-5]+360;
  end
  if strcmp(region,'sydney')==1
    lat         = -33.8;
    lon         = 150.0;
  end
  if strcmp(region,'dallas')==1
    lat         = 32.7;
    lon         = -96.8+360;
  end
  if strcmp(region,'seattle')==1
    lat         = 49.2;
    lon         = -121.4+360;
      var0      = 2;
      seas0      = 1;
  end
  if strcmp(region,'sahara')==1
    lat         = 22.5;
    lon         = 9.6;
      var0      = 2;
      seas0     = 3;
  end
  if strcmp(region,'denver')==1
    lat         = 39.9;
    lon         = -100.1+360;
  end
  if strcmp(region,'anchorage')==1
    lat         = 61.2;
    % lon         = -149.1+360;
    lon         = -147.4+360;
  end
  if strcmp(region,'north_america')==1
    lat         = [22.1 70.5];
    lon         = [-124.4+360 -53.7+360];
  end
  if strcmp(region,'kolkata')==1
    lat         = 22.6;
    lon         = 88.4;
  end
  if strcmp(region,'nao')==1
    lat         = 50; % dummy, actual NAO index is loaded later
    lon         = 50; % dummy, actual NAO index is loaded later
  end
  if strcmp(region,'southern_ocean')==1
    lat         = [-70 -60];
    lon         = [0 360];
      seas0      = 3;
      var0      = 1;
  end
  if strcmp(region,'sahel')==1
    lat         = [10 20];
    lon         = [0 40];
      var0      = 2;
      seas0     = 2;
  end
  if strcmp(region,'NH')==1
    lat         = [0 90];
    lon         = [0 360];
  end
  if strcmp(region,'SH')==1
    lat         = [-90 0];
    lon         = [0 360];
  end
  if strcmp(region,'alaska')==1
    lat         = [51.6 73.9];
    lon         = [191.2 255.9];
  end
  if strcmp(region,'zurich')==1
    lat         = 47.5;
    lon         = 8.5;
  end
  if strcmp(region,'arctic')==1
    lat         = [70 90];
    lon         = [0 360];
  end
  if strcmp(region,'AA')==1
    lat1        = [0 90];
    lat2        = [75 90];
    lon1        = [0 360];
    lon2        = [0 360];
  end
  % -- find grid cells to plot in models
  clear('tmp','tmp0')
  if strcmp(region,'global')==1
    ii      = 1:length(lati);
    jj      = 1:length(loni);
    weights = area(jj,ii);
    weights = weights/nansum(nansum(weights));
  else if strcmp(region,'Mediterranean')==1
    ii      = 1:length(lati);
    jj      = 1:length(loni);
    weights = area .* squeeze(srex_regions(:,:,12));
    weights = weights/nansum(nansum(weights));
  else if strcmp(region,'Central_Europe')==1
    ii      = 1:length(lati);
    jj      = 1:length(loni);
    weights = area .* squeeze(srex_regions(:,:,6));
    weights = weights/nansum(nansum(weights));
  else if strcmp(region,'Northern_Europe')==1
    ii      = 1:length(lati);
    jj      = 1:length(loni);
    weights = area .* squeeze(srex_regions(:,:,16));
    weights = weights/nansum(nansum(weights));
  else if strcmp(region,'NCA_nw')==1
    % -- NCA regions
    filein  = [pathin 'NCA_masks_g025.nc'];
    nca_region = ncread([filein],'NW');
    weights = area .* nca_region;
    weights = weights/nansum(nansum(weights));
  else if strcmp(region,'AA')==1
    for x = 1:length(lat1)
      tmp0    = find(abs(lati-lat1(x))==min(abs(lati-lat1(x))));
      tmp1(x)  = tmp0(1);
    end
    for x = 1:length(lat2)
      tmp0    = find(abs(lati-lat2(x))==min(abs(lati-lat2(x))));
      tmp2(x)  = tmp0(1);
    end
    ii1 = tmp1(1):tmp1(end);
    ii2 = tmp2(1):tmp2(end);
    for x = 1:length(lon1)
      tmp0    = find(abs(loni-lon1(x))==min(abs(loni-lon1(x))));
      tmp(x)  = tmp0(1);
    end
    jj1 = tmp(1):tmp(end);
    jj2 = jj1;
    % -- cut out area --
    weights1 = area(jj1,ii1);
    weights2 = area(jj2,ii2);
    weights1 = weights1/nansum(nansum(weights1));
    weights2 = weights2/nansum(nansum(weights2));
    clear('tmp1','tmp2')
  else
    for x = 1:length(lat)
      tmp0    = find(abs(lati-lat(x))==min(abs(lati-lat(x))));
      tmp(x)  = tmp0(1);
    end
    ii = tmp(1):tmp(end);
    for x = 1:length(lon)
      tmp0    = find(abs(loni-lon(x))==min(abs(loni-lon(x))));
      tmp(x)  = tmp0(1);
    end
    jj = tmp(1):tmp(end);
    if exist('lon2')==1
      clear('tmp')
      for x = 1:length(lon2)
        tmp0    = find(abs(loni-lon2(x))==min(abs(loni-lon2(x))));
        tmp(x)  = tmp0(1);
      end
      jj = [jj tmp(1):tmp(end)];
    end
    % -- cut out area --
    if strcmp(region,'nino34')==1 ||...
       strcmp(region,'southern_ocean')==1 ||...
       strcmp(region,'arctic')==1
      weights = area(jj,ii);
      weights = weights/nansum(nansum(weights));
    else
      weights = landarea(jj,ii);
      weights = weights/nansum(nansum(weights));
    end
  end
  end
  end
  end
  end
  end
  clear('tmp','tmp0')


  % -- variable --
  for v = [1] % <----------------------------------------------------------- EDIT!
    if r == 1
      vari = vars{v}
    else
      if exist('var0')==0
        var0 = v;
      end
      vari = vars{var0}
    end
    if strcmp(vari,'tas')==1
      units       = 'K';
      f           = 1;
      % units       = '\circF';
      % f           = 1.8;
      var_name    = 'Temperature'
      var_name2   = 'temperature';
    elseif strcmp(vari,'pr')==1
      units       = '%';
      f           = 86400;
      var_name    = 'Precipitation'
      var_name2   = 'precipitation';
    elseif strcmp(vari,'psl')==1
      units       = 'unitless';
      f           = 1;
      var_name    = 'Sea level pressure'
    end

    % -- load global mean observations
    obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202012_ts.nc'],'temperature'));
    obs_raw_am = am(obs_raw);
    if strcmp(vari,'tas')==1 && strcmp(region,'global')==1
      obs_name = 'BEST';
      time_obs = 1850:2022;
      % obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202012_ts.nc'],'temperature'));
      obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1_g025_ts.nc'],'temperature'));
      obs_raw_am = am(obs_raw);
      obs = rm(am(obs_raw),wl);
      obs = obs-nanmean(obs(refstart-1850+1:refende-1850+1));
    end

    % % -- load GMST (need always)
    % % time_tas_obs = 1850:2018;
    % % tas_obs_raw = ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-201812_ts.nc'],'TEMPERATURE');
    % time_tas_obs = 1850:2020;
    % tas_obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202012_ts.nc'],'temperature'));
    % tas_obs = rm(am(tas_obs_raw),wl);
    % ref_pi = nanmean(tas_obs(1850-1850+1:1899-1850+1)); % preindustrial reference
    % ref = nanmean(tas_obs(refstart-1850+1:refende-1850+1));
    % tas_obs_wtd = ref-ref_pi; % wtd = warming-to-date
    % tas_obs = tas_obs-nanmean(tas_obs(refstart-1850+1:refende-1850+1));




    % -- season --
    for s = [3] %1:length(seasons) % <------------------------------------------------------- EDIT!
      % seasons = {'DJF','JJA','annual'};

      % if r == 1
      %   seas = seasons{s};
      % else
      %   if exist('seas0')==0
      %     seas0 = s;
      %   end
      %   seas  = seasons{seas0};
      % end
      seas = seasons{s};

      % -- print to screen the current parameters
      '------------------------------------------------'
      region
      vari
      seas
      '------------------------------------------------'


      % -- load CMIP6 data --
      cmip6_ts_em  = NaN(length(models_cmip6),length(time));
      scen        = scen_cmip6;
      pathin_tmp  = [pathin 'cmip6-ng/' vari '/'];
      for sc = 1:length(scen)
        fid = fopen(['~/Dropbox/tmp.txt'], 'w');
        clear('ne')
        for m = 1:length(models_cmip6)
          ['scenario = ' scen{sc} ' / model = ' models_cmip6{m} ]
          files   = dir([pathin 'cmip6-ng/' vari '/' vari '_mon_' models_cmip6{m} '_' scen{sc} '_r*_g025.nc']);
          if sc == 4
            ne(m) = length(files);
          else
            ne(m) = 1;
          end
          for e = 1:ne(m)
            tmp0    = ncread([pathin_tmp files(e).name],[vari])*f; % yyy
            fprintf(fid, '%s\n', files(e).name);
            for i = 1:length(tmp0)
              if strcmp(region,'AA')==1
                tmp1(i)       = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
                tmp2(i)       = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
                tmp(i)     = tmp1(i)-tmp2(i);
              else
                tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
              end
            end
            tmp     = tmp((start-start_cmip6)*12+1:(ende-start_cmip6+1)*12);
            if strcmp(seas,'annual')==1
              tmp    = am(tmp);
            else
              tmp    = seasmean(tmp,seas);
            end
            eval(['cmip6_' scen{sc} '_ts_raw(m,e,:) = tmp;'])
            ref = nanmean(tmp(refstart-start+1:refende-start+1));
            if strcmp(vari,'pr')==1
              eval(['cmip6_' scen{sc} '_ts(m,:)   = ((tmp-ref)/ref)*100;'])
            else
              eval(['cmip6_' scen{sc} '_ts(m,:)   = tmp-ref;'])
            end
          end
          % -- H&S style calculation --
          eval(['tmp1 = squeeze(cmip6_' scen{sc} '_ts_raw(m,1,:));'])
          idx = ~isnan(tmp1);
          fit = NaN(size(tmp1));
          fit(idx) = polyval(polyfit(time(idx),tmp1(idx)',4),time(idx)); % polynomial fit
          idx = ~isnan(fit);
          ref = nanmean(fit(refstart-start+1:refende-start+1));
          if strcmp(vari,'pr')==1
            eval(['cmip6_' scen{sc} '_ts_em(m,:)   = ((fit-ref)/ref)*100;'])
          else
            eval(['cmip6_' scen{sc} '_ts_em(m,:)   = fit-ref;'])
          end
          tmp1 = tmp1';
          residual  = tmp1-fit' + nanmean(fit);
          if strcmp(vari,'pr')==1
            ref       = nanmean(residual);
            residual  = (residual/ref)*100;
          end
          tmp2 = rm(residual,wl);
          eval(['cmip6_' scen{sc} '_ts_residual(m,:)  = tmp2;'])
          eval(['cmip6_' scen{sc} '_ts_var1(m)        = nanvar(tmp2);'])
          eval(['cmip6_' scen{sc} '_ts_var2(m)        = nanvar(tmp2(1:refende-start+1));'])
          eval(['cmip6_' scen{sc} '_ts_var3(m)        = nanvar(tmp2(refende-start+1+1:end));'])
        end
        fclose(fid);
      end


      if load_dple == 1
        tic
        % -- CESM1-DPLE
        'load CESM1-DPLE...'
        % -- DPLE anom (de-drifted, which is effectively also de-seasonalized)
        filein  = [pathin '/DPLE/Amon/' vari '/' vari '_Amon_CESM1-DPLE_g025_anom.nc'];
        tmp0    = ncread([filein],[vari])*f;
        tmp     = NaN(40,122,64);
        for e = 1:40
          e
          for lt = 1:122
            for y = 1:64
              tmp(e,lt,y)     = nansum(nansum(squeeze(tmp0(jj,ii,e,lt,y)).*weights));
            end
          end
        end
        clear('tmp0','tmp2','tmp3')
        % tmp = tmp(:,3:end,:); % cut Nov+Dec to start in Jan
        % tmp = tmp((start-1954)*12+1:(ende-1954+1)*12);
        for e = 1:40
          for y = 1:64
            if strcmp(seas,'annual')==1
              tmp2(e,:,y)    = am(squeeze(tmp(e,3:end,y)));
            elseif strcmp(seas,'DJF')==1
              for lt = 1:floor(122/12)
                tmp99 = squeeze(tmp(e,(lt-1)*12+2:(lt-1)*12+4,y));
                tmp2(e,lt,y) = tmp99(1)*(31/90) + tmp99(2)*(31/90) + tmp99(3)*(28/90);
              end
            else
              tmp2(e,:,y)    = seasmean(squeeze(tmp(e,:,y)),seas);
            end
            if wl > 1
              tmp3(e,:,y) = rm(tmp2(e,:,y),wl);
            end
          end
        end
        dple      = tmp2 + nanmean(obs(1955-time_obs(1)+1:2017-time_obs(1))); % correct for different reference periods
        save(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/dple_' region '_' seas '.mat'], 'dple');
        if wl > 1
          dple_rm = tmp3;
          save(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/dple_' region '_' seas '_rm.mat'], 'dple_rm');
        end
        % -- DPLE raw ----
        filein  = [pathin '/DPLE/Amon/' vari '/' vari '_Amon_CESM1-DPLE_g025.nc'];
        tmp0    = ncread([filein],[vari])*f;
        tmp     = NaN(40,122,64);
        for e = 1:40
          e
          for lt = 1:122
            for y = 1:64
              tmp(e,lt,y)     = nansum(nansum(squeeze(tmp0(jj,ii,e,lt,y)).*weights));
            end
          end
        end
        clear('tmp0','tmp2','tmp3')
        for e = 1:40
          for y = 1:64
            if strcmp(seas,'annual')==1
              tmp2(e,:,y)    = am(squeeze(tmp(e,3:end,y)));
            elseif strcmp(seas,'DJF')==1
              for lt = 1:floor(122/12)
                tmp99 = squeeze(tmp(e,(lt-1)*12+2:(lt-1)*12+4,y));
                tmp2(e,lt,y) = tmp99(1)*(31/90) + tmp99(2)*(31/90) + tmp99(3)*(28/90);
              end
            else
              tmp2(e,:,y)    = seasmean(squeeze(tmp(e,:,y)),seas);
            end
            if wl > 1
              tmp3(e,:,y) = rm(tmp2(e,:,y),wl);
            end
          end
        end
        dple_raw = tmp2;
        save(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/dple_raw_' region '_' seas '.mat'], 'dple_raw');
        if wl > 1
          dple_rm = tmp3;
          save(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/dple_' region '_' seas '_rm.mat'], 'dple_raw_rm');
        end
        % ref = nanmean(tmp(refstart-start+1:refende-start+1));
        % if strcmp(vari,'pr')==1
        %   tmp   = ((tmp-ref)/ref)*100;
        % else
        %   tmp   = tmp-ref;
        % end
        % dple_ensmean(sc,:) = tmp;
        % idx = ~isnan(tmp);
        % dple_ensmean2(sc,:) = polyval(polyfit(time(idx),tmp(idx),4),time);
        toc
      else
        % dple_raw  = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/dple_raw_' region '_' seas '.mat']).dple_raw;
        % dple      = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/dple_' region '_' seas '.mat']).dple;
      end


    end % -- end of seasons loop
  end % -- end of vars loop
end % -- end of regions loop




      %% -- UNCERTAINTIES --
      % -- normal HS09 uncertainties --
      % -- CMIP6 scen uncertainty with 4 scenarios:
      scen_u_cmip6      = nanvar([nanmean(cmip6_ssp585_ts_em,1); nanmean(cmip6_ssp370_ts_em,1); nanmean(cmip6_ssp245_ts_em,1); nanmean(cmip6_ssp126_ts_em,1)]); % using 4 scenarios

      % -- observed warming since preindustrial (today = reference period)
      if strcmp(vari,'tas')==1
        if strcmp(region,'global')==1
          offset = nanmean(obs_raw_am(refstart-time_obs(1)+1:refende-time_obs(1))) - nanmean(obs_raw_am(1880-time_obs(1)+1:1910-time_obs(1)))
        else
          offset = nanmean(obs(refstart-time_obs(1)+1:refende-time_obs(1))) - nanmean(obs(1880-time_obs(1)+1:1910-time_obs(1)))
        end
      else
        offset = 0;
      end

      [a,i] = min(scen_u_cmip6);
      scen_u_cmip6(1:i) = 0;
      scen_u_cmip6_3scen(1:i) = 0;
      scen_u_cmip6_brekke(1:i) = 0;

      % -- what HS09 did, namely calculate the mean across model uncertainties from the different scenarios):
      model_u_cmip6_hs  = nanmean([nanvar(cmip6_ssp585_ts_em,1); nanvar(cmip6_ssp370_ts_em,1); nanvar(cmip6_ssp245_ts_em,1); nanvar(cmip6_ssp126_ts_em,1)]);
      model_u_cmip6     = model_u_cmip6_hs;
      int_u_cmip6_mean  = nanmean(cmip6_ssp585_ts_var1);
      int_u_cmip6_max   = max(cmip6_ssp585_ts_var1);
      int_u_cmip6_min   = min(cmip6_ssp585_ts_var1);
      total_u_cmip6_mean        = scen_u_cmip6 + model_u_cmip6 + int_u_cmip6_mean;
      total_u_cmip6_mean_hs     = scen_u_cmip6 + model_u_cmip6_hs + int_u_cmip6_mean;

      % -- fractional uncertainties --
      scen_u_frac_cmip6     = (scen_u_cmip6./total_u_cmip6_mean)*100;
      scen_u_frac_cmip6_hs  = (scen_u_cmip6./total_u_cmip6_mean_hs)*100;
      model_u_frac_cmip6    = (model_u_cmip6./total_u_cmip6_mean)*100;
      model_u_frac_cmip6_hs = (model_u_cmip6_hs./total_u_cmip6_mean_hs)*100;
      int_u_frac_cmip6_mean     = (int_u_cmip6_mean./total_u_cmip6_mean)*100;


      % % -- load global mean temperature change (forced response, NOTE: added a warming-to-date "tas_obs_wtd" estimate from observations):
      % % offset = tas_obs_wtd;
      % % offset = 0;
      % tic
      % tasg_le_ts_em            = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_le_' num2str(wl) 'yr.txt'])' + offset;
      % tasg_cmip5_rcp26_ts_em   = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_cmip5_rcp26_em_' num2str(wl) 'yr.txt'])' + offset;
      % tasg_cmip5_rcp45_ts_em   = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_cmip5_rcp45_em_' num2str(wl) 'yr.txt'])' + offset;
      % tasg_cmip5_rcp85_ts_em   = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_cmip5_rcp85_em_' num2str(wl) 'yr.txt'])' + offset;
      % tasg_cmip6_ssp126_ts_em  = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_cmip6_ssp126_em_' num2str(wl) 'yr.txt'])' + offset;
      % tasg_cmip6_ssp245_ts_em  = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_cmip6_ssp245_em_' num2str(wl) 'yr.txt'])' + offset;
      % tasg_cmip6_ssp370_ts_em  = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_cmip6_ssp370_em_' num2str(wl) 'yr.txt'])' + offset;
      % tasg_cmip6_ssp585_ts_em  = load(['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/time_series/tas_global_cmip6_ssp585_em_' num2str(wl) 'yr.txt'])' + offset;
      % % -- sort variable according to tasg:
      % tmp   = [tasg_le_ts_em(:); tasg_cmip5_rcp85_ts_em(:); tasg_cmip6_ssp585_ts_em(:)] + offset;
      % range = [min(tmp):(max(tmp)-min(tmp))/100:max(tmp)] + offset;
      % clear('le_ts_em_vs_tasg','le_var_vs_tasg','cmip5_rcp85_ts_em_vs_tasg','cmip6_ssp585_ts_em_vs_tasg')
      % le_ts_em_vs_tasg            = NaN(length(models),length(range));
      % le_var_vs_tasg              = NaN(length(models),length(range));
      % cmip5_rcp85_ts_em_vs_tasg   = NaN(length(models_cmip5),length(range));
      % cmip6_ssp585_ts_em_vs_tasg  = NaN(length(models_cmip6),length(range));
      % for i = 1:length(range)
      %   for m = 1:length(models)
      %     [null,idx] = min(abs(tasg_le_ts_em(:,m)-range(i)));
      %     le_ts_em_vs_tasg(m,i) = le_ts_em(m,idx);
      %     le_var_vs_tasg(m,i) = nanvar(le_ts_tmp(m,:,idx));
      %     % le_var_vs_tasg(m,i) = int_u_le_mean(idx);
      %   end
      %   for s = 1:length(scen_cmip5)
      %     for m = 1:length(models_cmip5)
      %       eval(['tmp = islocalmin((abs(tasg_cmip5_' scen_cmip5{s} '_ts_em(:,m)-range(i))));']);
      %       tmp = find(tmp==1);
      %       if isempty(tmp)==1
      %         eval(['[null,idx] = min(abs(tasg_cmip5_' scen_cmip5{s} '_ts_em(:,m)-range(i)));']);
      %       else
      %         idx = tmp(1);
      %       end
      %       eval(['cmip5_' scen_cmip5{s} '_ts_em_vs_tasg(m,i) = cmip5_' scen_cmip5{s} '_ts_em(m,idx);']);
      %     end
      %   end
      %   for s = 1:length(scen_cmip6)
      %     for m = 1:length(models_cmip6)
      %       eval(['[null,idx] = min(abs(tasg_cmip6_' scen_cmip6{s} '_ts_em(:,m)-range(i)));']);
      %       eval(['cmip6_' scen_cmip6{s} '_ts_em_vs_tasg(m,i) = cmip6_' scen_cmip6{s} '_ts_em(m,idx);']);
      %     end
      %   end
      % end

      % model_u_cmip6_ssp126_tasg = nanvar(cmip6_ssp126_ts_em_vs_tasg);
      % model_u_cmip6_ssp245_tasg = nanvar(cmip6_ssp245_ts_em_vs_tasg);
      % model_u_cmip6_ssp370_tasg = nanvar(cmip6_ssp370_ts_em_vs_tasg);
      % model_u_cmip6_ssp585_tasg = nanvar(cmip6_ssp585_ts_em_vs_tasg);
      % int_u_cmip6_tasg_mean     = int_u_cmip6_mean;
      % total_u_cmip6_ssp126_tasg_mean   = model_u_cmip6_ssp126_tasg + int_u_cmip6_tasg_mean;
      % total_u_cmip6_ssp245_tasg_mean   = model_u_cmip6_ssp245_tasg + int_u_cmip6_tasg_mean;
      % total_u_cmip6_ssp370_tasg_mean   = model_u_cmip6_ssp370_tasg + int_u_cmip6_tasg_mean;
      % total_u_cmip6_ssp585_tasg_mean   = model_u_cmip6_ssp585_tasg + int_u_cmip6_tasg_mean;
      %
      % model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
      % model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
      % model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
      % model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
      % int_u_frac_cmip6_ssp126_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp126_tasg_mean)*100;
      % int_u_frac_cmip6_ssp245_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp245_tasg_mean)*100;
      % int_u_frac_cmip6_ssp370_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp370_tasg_mean)*100;
      % int_u_frac_cmip6_ssp585_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp585_tasg_mean)*100;





      % --- PLOTTING ---------------------------------------------------------------
      close all

      cols = [255 0 0;
              255 160 16;
              255 224 32;
              0 192 0;
              80 208 255;
              0 32 255;
              160 32 255]/255;

      hs09_cols = ...
      [53 74 161;...
       255 110 4;...
       0 127 60]/255;
      hs09_cols_light = ...
      [164 180 245;...
       252 210 179;...
       172 232 200]/255;

      rcp_cols = [217,37,42;...
                 232,126,63;...
                 133,177,212;...
                 51,74,141]/255;
      rcp_cols_light = ...
      [235,144,115;...
      243,183,136;...
      189,211,227;...
      138,141,185]/255;
      % rcp_cols_light2 = (1-rcp_cols)*.25+rcp_cols;

      if wl > 1
        xlim  = [2015 2099-round(wl/2)];
        xlim0 = [1950+round(wl/2) 2099-round(wl/2)];
      else
        xlim  = [refende+1 2099-round(wl/2)];
        xlim0 = [1950+round(wl/2) 2099-round(wl/2)];
      end
      if strcmp(vari,'tas')==1
        yincr = 1;
        if refstart == 1995 || refstart == 2001
          ylim = [-1.25 5.5]*f;
        else
          ylim = [0 5.5];
        end
      else % vari = pr
        yincr = 2;
        ylim = [-2.9 10];
        if strcmp(region,'Mediterranean')==1
          ylim = [-10 10];
        end
        if strcmp(region,'UCRB')==1
          ylim = [-10 20];
        end
      end

      tl = 0.015; % tick length
      sh = 0.05; % horizontal space between panels
      sv = 0.04; % vertical space between panels




      % ----------------------------------
      if fig3 == 1;
      % -- plot --------
      close all

        ylim = [-1.8 6.5];
        if wl == 1
          xlim = [xlim(1)-10 xlim(2)];
        end

        figure1 = figure;
        set(figure1, 'units', 'centimeters', 'pos', [10 10 26 6]) % for 3-row final paper figure

        % -- CMIP6:
        subaxis(1,3,1, 'sh', sh, 'sv', sv)
        hold on
        title(['CMIP6 projections rel. to ' num2str(refstart) '-' num2str(refende)],'Interpreter','none')
        plot(time,cmip6_ssp585_ts','Color',rcp_cols_light(1,:))
        plot(time,cmip6_ssp370_ts','Color',rcp_cols_light(2,:))
        plot(time,cmip6_ssp245_ts','Color',rcp_cols_light(3,:))
        plot(time,cmip6_ssp126_ts','Color',rcp_cols_light(4,:))
        % plot(time,cmip6_ssp585_ts_em','Color',rcp_cols_light(1,:))
        % plot(time,cmip6_ssp370_ts_em','Color',rcp_cols_light(2,:))
        % plot(time,cmip6_ssp245_ts_em','Color',rcp_cols_light(3,:))
        % plot(time,cmip6_ssp126_ts_em','Color',rcp_cols_light(4,:))
        h4 = plot(time,nanmean(cmip6_ssp585_ts_em),'Color',rcp_cols(1,:),'LineWidth',3)
        h3 = plot(time,nanmean(cmip6_ssp370_ts_em),'Color',rcp_cols(2,:),'LineWidth',3)
        h2 = plot(time,nanmean(cmip6_ssp245_ts_em),'Color',rcp_cols(3,:),'LineWidth',3)
        h1 = plot(time,nanmean(cmip6_ssp126_ts_em),'Color',rcp_cols(4,:),'LineWidth',3)
        % if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
        %   h0 = plot(time_obs,obs,'k','LineWidth',2)
        % end
        if wl == 1 && strcmp(region,'anchorage')==1 && strcmp(seas,'DJF')==1
          set(gca,'YLim',[-9 15],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
        else
          set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
        end
        ylim = get(gca,'YLim');
        hline(0,'k')
        box on
        % legend([h3 h2 h1],'RCP8.5','RCP4.5','RCP2.6','Location','NorthWest','FontSize',9)
        legend([h4 h3 h2 h1],'SSP5-8.5','SSP3-7.0','SSP2-4.5','SSP1-2.6','Location','NorthWest','FontSize',9)
        legend boxoff
        xlabel('Time (Year)')
        ylabel(['Temperature change (' units ')'])
        if r==1
          text(xlim(1)+.8*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(a)')
        else
          text(xlim(1)+.8*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(d)')
        end

        subaxis(1,3,2, 'sh', sh, 'sv', sv)
        hold on
        % title(['Uncertainty for' char(10) region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units '^2)'],'Interpreter','none')%,'FontSize',10)
        title(['Sources of uncertainty'],'Interpreter','none')%,'FontSize',10)
        % -- variance
        h1 = plot(time,model_u_cmip6*0+int_u_cmip6_mean,'Color',[255 110 4]/255,'LineWidth',2);
        % h1 = plot(time,rm(int_u_le_mean,10),'Color',[255 110 4]/255,'LineWidth',2);
        h2 = plot(time,model_u_cmip6,'Color',[53 74 161]/255,'LineWidth',2);
        h3 = plot(time,scen_u_cmip6,'Color',[0 127 60]/255,'LineWidth',2);
        h0 = plot(time,total_u_cmip6_mean,'Color','k','LineWidth',2);
        % set(gca,'YLim',[0 ylim(end)],'YTick',[0:yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
        if wl == 1 && strcmp(region,'anchorage')==1 && strcmp(seas,'DJF')==1
          set(gca,'YLim',[0 15],'XLim',xlim,'TickLength',[tl tl],'Layer','top')
        else
          set(gca,'YLim',[0 3],'YTick',[0:yincr:round(ylim(end))],'XLim',xlim,'TickLength',[tl tl],'Layer','top')
        end
        ylim = get(gca,'YLim');
        box on
        legend([h0(1) h1(1) h2(1) h3(1)],'Total','Internal variability','Response','Scenario','Location','NorthWest','FontSize',9,'Interpreter','none')
        legend boxoff
        xlabel('Time (Year)')
        ylabel(['Variance (' units '^2)'])
        if r==1
          % text(2085,0+abs(.08*(3-0)),'(b)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(b)')
        else
          % text(2085,0+abs(.08*(3-0)),'(e)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(e)')
        end

        subaxis(1,3,3, 'sh', sh, 'sv', sv)
        hold on
        title({'Contribution to total uncertainty'})%,'FontSize',10)
        idx       = find(~isnan(total_u_cmip6_mean));
        x         = time(idx);
        y0        = x*0;
        ym_cmip6  = model_u_frac_cmip6(idx); % model u
        yms_cmip6 = model_u_frac_cmip6(idx) + scen_u_frac_cmip6(idx); % model u plus scen u
        ytop      = time(idx)*0+100;
        h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip6)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
        h1 = patch([x fliplr(x)],[yms_cmip6 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
        h3 = patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
        patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
        hold on
        set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
        ylim = get(gca,'YLim');
        if r==1
          legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario','Response','Location','best','FontSize',9)
        end
        ylabel({'Fraction (%)'})
        xlabel('Time (Year)')
        box on
        if r==1
          % text(2085,0+abs(.08*(100)),'(c)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(c)')
        else
          % text(2085,0+abs(.08*(100)),'(f)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(f)')
        end

        tightfig

        % return
        % set(gcf,'PaperPositionMode','auto');
        % set(gcf,'renderer','Painters')
        % fileo = ['~/Dropbox/publication/lehner23_intvar_perspective/fig/hawkins_plots_1x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_cmip6'];
        % print('-r300','-loose', '-depsc', ['' fileo '.eps'])
        % save2pdf(['' fileo '.pdf'])
        % saveas(gcf,fileo,'jpg')
        % return

      end



      % ----------------------------------
      if fig5 == 1;
      % -- plot --------
      close all

        % xlim = [1950 2100];
        xlim = [2000 2030];

        figure1 = figure;
        % set(figure1, 'units', 'centimeters', 'pos', [10 10 26 7]) % for 2-row final paper figure
        set(figure1, 'units', 'centimeters', 'pos', [10 10 26 6]) % for 3-row final paper figure

        subplot(1,3,1)
        hold on
        ny = 5; % 6% how many years before last forecast init (set to 0 if last forecast init is desired, which is from Nov 2017)
        % plot(time,squeeze(le_raw_ts(1,:,:))-nanmean(nanmean(le_raw_ts(1,:,1955-start+1:2017-start),2)),'Color',[.7 .7 .7])
        plot(time,squeeze(le_ts_anom(1,1:35,:)),'Color',[.7 .7 .7])
        plot(2018-ny:2027-ny,squeeze(dple(:,:,end-ny)),'Color',[.5 .5 1])
        h1 = plot(time,nanmean(squeeze(le_ts_anom(1,1:35,:))),'k','LineWidth',1)
        h2 = plot(2018-ny:2027-ny,squeeze(nanmean(dple(:,:,end-ny),1)),'b','LineWidth',2)
        h3 = plot(time_obs,obs,'k.-','LineWidth',2,'MarkerSize',12)
        ylim = get(gca,'YLim');
        if strcmp(regions{r},'global')==1 && strcmp(seas,'annual')==1
          ylim = [-.5 1];
        end
        if strcmp(regions{r},'north_america')==1 && strcmp(seas,'DJF')==1
          ylim = [-4 5];
        end
        if strcmp(regions{r},'anchorage')==1 && strcmp(seas,'DJF')==1
          ylim = [-8 10];
        end
        set(gca,'Layer','top','XLim',xlim,'YLim',ylim)
        hline(0,'k')
        box on
        ylabel('Temperature anomaly (\circC)')
        xlabel('Time (Year)')
        title([regions{r} ' ' seas ' temperature'])
        legend([h1 h2 h3],'CESM1-LE','CESM1-DPLE','Observations','location','northwest')
        legend boxoff
        if strcmp(region,'global')==1
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(a)')
        else
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(d)')
        end


        subplot(1,3,2)
        hold on
        x = 1:10;
        xlim2 = [x(1) x(end)];
        clear('tmp1','tmp2')
        % tmp1 = nanstd(squeeze(le_ts_anom(1,1:35,1955-1920+1:2017-1920)),1);
        tmp1 = nanstd(squeeze(le_raw_ts(1,1:35,1955-1920+1:2017-1920)),1);
        jbfill(x,x*0+max(tmp1),x*0+min(tmp1),[.9 .9 .9],'none')
        jbfill(x,x*0+prctile(tmp1,5),x*0+prctile(tmp1,95),[.7 .7 .7],'none')
        hold on
        h1 = plot(x,x*0+nanmean(tmp1),'k','LineWidth',1)
        for i = 1:64
          % tmp2(i,:) = nanstd(squeeze(dple(:,:,i)));
          tmp2(i,:) = nanstd(squeeze(dple_raw(1:35,:,i)));
          % plot(x,tmp2(i,:),'Color',[40, 86, 237]/255)
        end
        jbfill(x,x*0+max(tmp2),x*0+min(tmp2),[.7 .7 1],'none',1,.5)
        jbfill(x,x*0+prctile(tmp2,5),x*0+prctile(tmp2,95),[.5 .5 1],'none',1,.5)
        hold on
        h2 = plot(x,nanmean(tmp2),'b','LineWidth',2)
        ratio = nanmean(tmp2)./nanmean(tmp1);
        ratio(ratio>1) = 1;
        set(gca,'Layer','top','XLim',xlim2,'XTick',[xlim2(1):1:xlim2(2)])
        xtickangle(0)
        ylim2 = get(gca,'YLim');
        % set(gca,'YTick',[ceil(ylim2(1)):.02:round(ylim2(end))])
        box on
        ylabel('Ensemble spread, \sigma (\circC)')
        xlabel('Lead time (Years)')
        title('Influence of initialization')
        legend([h1 h2],'CESM1-LE','CESM1-DPLE','location','south')
        legend boxoff
        if strcmp(region,'global')==1
          text(xlim2(1)+.88*(xlim2(2)-xlim2(1)),ylim2(1)+abs(.1*(ylim2(2)-ylim2(1))),'(b)')
        else
          text(xlim2(1)+.88*(xlim2(2)-xlim2(1)),ylim2(1)+abs(.1*(ylim2(2)-ylim2(1))),'(e)')
        end


        subplot(1,3,3)
        hold on
        if constrain == 0
          tmp = nanmean([nanmean(cmip6_ssp585_ts_em,1); nanmean(cmip6_ssp370_ts_em,1); nanmean(cmip6_ssp245_ts_em,1); nanmean(cmip6_ssp126_ts_em,1)]);
          idx = ~isnan(tmp);
          i = sqrt(int_u_cmip6_mean);
          % m = sqrt(model_u_cmip6_hs(idx));
          m = sqrt(model_u_cmip6(idx));
          s = sqrt(scen_u_cmip6(idx));
        else
          tmp = nanmean([nanmean(weights_cmip6.*cmip6_ssp585_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp370_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp245_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp126_ts_em,1)]);
          idx = ~isnan(tmp);
          i = sqrt(int_u_cmip6_w_mean);
          m = sqrt(model_u_cmip6_w(idx));
          s = sqrt(scen_u_cmip6_w(idx));
        end
        % -- build in constrained int var from DPLE
        int = s+m+i;
        [intmin,idx2] = min(int);
        if strcmp(seas,'DJF')==1
          idx2 = idx2+2;
        else
          idx2 = idx2+1;
        end
        int_constr = int;
        int_constr(idx2:idx2+9) = int(idx2:idx2+9).*(ratio);
        int_constr(idx2-1) = 0;
        % --
        h3 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(s+m+i) fliplr(tmp(idx)+(s+m+i))],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
        h2 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(s+m) fliplr(tmp(idx)+(s+m))],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
        h1 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-s fliplr(tmp(idx)+s)],hs09_cols(3,:),'LineWidth',.1,'Edgecolor','none')
        h4 = plot(time(idx2-1:end),tmp(idx2-1:end)+(int_constr(idx2-1:end)),'k:','LineWidth',2)
        plot(time(idx2-1:end),tmp(idx2-1:end)-(int_constr(idx2-1:end)),'k:','LineWidth',2)
        hold on
        plot(time_obs,obs,'k.-','LineWidth',2,'MarkerSize',12)
        % set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim,'TickLength',[tl tl],'Layer','top')
        set(gca,'YLim',ylim,'XLim',xlim,'Layer','top')
        hline(0,'k')
        box on
        ylabel('Temperature anomaly (\circC)')
        xlabel('Time (Year)')
        title('Source of uncertainty')
        % legend([h3(1) h4(1) h2(1) h1(1)],'Internal variability', ['Internal variability' char(10) 'constrained'],'Response','Scenario','Location','NorthWest','Interpreter','none')
        legend([h3(1) h4(1) h2(1) h1(1)],'Internal variability', ['Int. var. constrained'],'Response','Scenario','Location','NorthWest','Interpreter','none')
        legend boxoff
        if strcmp(region,'global')==1
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(c)')
        else
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(f)')
        end

        tightfig

        return
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'renderer','Painters')
        % fileo = ['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/fig/hawkins_plots/' vari '/hawkins_plots_1x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_dple'];
        fileo = ['~/Dropbox/publication/lehner23_intvar_perspective/fig/hawkins_plots_1x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_dple'];
        print('-r300','-loose', '-depsc', ['' fileo '.eps'])
        save2pdf(['' fileo '.pdf'])
        % saveas(gcf,fileo,'jpg')
      end



      % ----------------------------------
      if fig3 == 1;
      % -- plot --------
      close all

        ylim = [-1.8 6.5];
        if wl == 1
          xlim = [xlim(1)-10 xlim(2)];
        end

        figure1 = figure;
        % set(figure1, 'units', 'centimeters', 'pos', [10 10 26 7]) % for 2-row final paper figure
        set(figure1, 'units', 'centimeters', 'pos', [10 10 26 6]) % for 3-row final paper figure

        % -- CMIP6:
        subaxis(1,3,1, 'sh', sh, 'sv', sv)
        hold on
        % title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means' char(10) 'relative to ' num2str(refstart) '-' num2str(refende)],'Interpreter','none')%,'FontSize',10)
        title(['CMIP6 projections rel. to ' num2str(refstart) '-' num2str(refende)],'Interpreter','none')
        plot(time,cmip6_ssp585_ts','Color',rcp_cols_light(1,:))
        plot(time,cmip6_ssp370_ts','Color',rcp_cols_light(2,:))
        plot(time,cmip6_ssp245_ts','Color',rcp_cols_light(3,:))
        plot(time,cmip6_ssp126_ts','Color',rcp_cols_light(4,:))
        % plot(time,cmip6_ssp585_ts_em','Color',rcp_cols_light(1,:))
        % plot(time,cmip6_ssp370_ts_em','Color',rcp_cols_light(2,:))
        % plot(time,cmip6_ssp245_ts_em','Color',rcp_cols_light(3,:))
        % plot(time,cmip6_ssp126_ts_em','Color',rcp_cols_light(4,:))
        h4 = plot(time,nanmean(cmip6_ssp585_ts_em),'Color',rcp_cols(1,:),'LineWidth',3)
        h3 = plot(time,nanmean(cmip6_ssp370_ts_em),'Color',rcp_cols(2,:),'LineWidth',3)
        h2 = plot(time,nanmean(cmip6_ssp245_ts_em),'Color',rcp_cols(3,:),'LineWidth',3)
        h1 = plot(time,nanmean(cmip6_ssp126_ts_em),'Color',rcp_cols(4,:),'LineWidth',3)
        % if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
        %   h0 = plot(time_obs,obs,'k','LineWidth',2)
        % end
        if wl == 1 && strcmp(region,'anchorage')==1 && strcmp(seas,'DJF')==1
          set(gca,'YLim',[-9 15],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
        else
          set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
        end
        ylim = get(gca,'YLim');
        hline(0,'k')
        box on
        % legend([h3 h2 h1],'RCP8.5','RCP4.5','RCP2.6','Location','NorthWest','FontSize',9)
        legend([h4 h3 h2 h1],'SSP5-8.5','SSP3-7.0','SSP2-4.5','SSP1-2.6','Location','NorthWest','FontSize',9)
        legend boxoff
        xlabel('Time (Year)')
        ylabel(['Temperature change (' units ')'])
        if r==1
          text(xlim(1)+.8*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(a)')
        else
          text(xlim(1)+.8*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(d)')
        end

        subaxis(1,3,2, 'sh', sh, 'sv', sv)
        hold on
        % title(['Uncertainty for' char(10) region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units '^2)'],'Interpreter','none')%,'FontSize',10)
        title(['Sources of uncertainty'],'Interpreter','none')%,'FontSize',10)
        % -- variance
        h1 = plot(time,model_u_cmip6*0+int_u_cmip6_mean,'Color',[255 110 4]/255,'LineWidth',2);
        % h1 = plot(time,rm(int_u_le_mean,10),'Color',[255 110 4]/255,'LineWidth',2);
        h2 = plot(time,model_u_cmip6,'Color',[53 74 161]/255,'LineWidth',2);
        h3 = plot(time,scen_u_cmip6,'Color',[0 127 60]/255,'LineWidth',2);
        h0 = plot(time,total_u_cmip6_mean,'Color','k','LineWidth',2);
        % set(gca,'YLim',[0 ylim(end)],'YTick',[0:yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
        if wl == 1 && strcmp(region,'anchorage')==1 && strcmp(seas,'DJF')==1
          set(gca,'YLim',[0 15],'XLim',xlim,'TickLength',[tl tl],'Layer','top')
        else
          set(gca,'YLim',[0 3],'YTick',[0:yincr:round(ylim(end))],'XLim',xlim,'TickLength',[tl tl],'Layer','top')
        end
        ylim = get(gca,'YLim');
        box on
        legend([h0(1) h1(1) h2(1) h3(1)],'Total','Internal variability','Response','Scenario','Location','NorthWest','FontSize',9,'Interpreter','none')
        legend boxoff
        xlabel('Time (Year)')
        ylabel(['Variance (' units '^2)'])
        if r==1
          % text(2085,0+abs(.08*(3-0)),'(b)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(b)')
        else
          % text(2085,0+abs(.08*(3-0)),'(e)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(e)')
        end

        subaxis(1,3,3, 'sh', sh, 'sv', sv)
        hold on
        % title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
        title({'Contribution to total uncertainty'})%,'FontSize',10)
        idx       = find(~isnan(total_u_cmip6_mean));
        x         = time(idx);
        y0        = x*0;
        % ym_cmip6  = model_u_frac_cmip6_hs(idx); % model u
        % yms_cmip6 = model_u_frac_cmip6_hs(idx) + scen_u_frac_cmip6_hs(idx); % model u plus scen u
        if constrain == 0
          ym_cmip6  = model_u_frac_cmip6(idx); % model u
          yms_cmip6 = model_u_frac_cmip6(idx) + scen_u_frac_cmip6(idx); % model u plus scen u
        else
          ym_cmip6  = model_u_frac_cmip6_w(idx); % model u
          yms_cmip6 = model_u_frac_cmip6_w(idx) + scen_u_frac_cmip6_w(idx); % model u plus scen u
        end
        ytop      = time(idx)*0+100;
        h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip6)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
        h1 = patch([x fliplr(x)],[yms_cmip6 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
        h3 = patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
        patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
        hold on
        set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
        ylim = get(gca,'YLim');
        if r==1
          legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario','Response','Location','best','FontSize',9)
        end
        ylabel({'Fraction (%)'})
        xlabel('Time (Year)')
        box on
        if r==1
          % text(2085,0+abs(.08*(100)),'(c)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(c)')
        else
          % text(2085,0+abs(.08*(100)),'(f)')
          text(xlim(1)+.88*(xlim(2)-xlim(1)),ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(f)')
        end

        tightfig

        return
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'renderer','Painters')
        % fileo = ['~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/fig/hawkins_plots/' vari '/hawkins_plots_1x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_cmip6'];
        fileo = ['~/Dropbox/publication/lehner23_intvar_perspective/fig/hawkins_plots_1x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_cmip6'];
        print('-r300','-loose', '-depsc', ['' fileo '.eps'])
        save2pdf(['' fileo '.pdf'])
        saveas(gcf,fileo,'jpg')
        % return

      end


%     end % -- end of seasons loop
%   end % -- end of vars loop
% end % -- end of regions loop



% -- functions -------------
function model_resp = mr(m,x,f,sens)
  % -- model, time series, forcing, sensitivity
  model_resp = x+(f*sens(m)); %.^.8;
  % model_resp = x.*(f);
end
