load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "~asphilli/ncl/asphilli_util.ncl"

; ---------------------------
; lehner23erc_fig4.ncl
; ---------------------------
; Code to produce Fig. 4 in Lehner and Deser (2023)
;
; Citation:
; Lehner, F., C. Deser (2023):
; Origin, importance, and predictive limits of internal climate variability
; Environmental Research: Climate, DOI: https://doi.org/10.1088/2752-5295/accf30
;
; Notes on code:
; - Only does calculations and plotting, not pre-processing
; - Requires concatenated monthly mean data from:
;   - CESM1-LE (b.e11.B20TRC5CNBDRD.f09_g16.*.cam.h0.*.nc and b.e11.BRCP85C5CNBDRD.f09_g16.*.cam.h0.*.nc)
;   - Dynamical component and residual component from CESM1-LE (the dynamical adjustment code can be found here: https://github.com/terrayl/Dynamico)
; ---------------------------

begin

  ; -- parameters --

  recalc1 = "yes" ; recalc CESM trends?
  recalc2 = "yes" ; recalc obs trends?

  ; -- season
  seas  = "DJF"; wateryear annual DJF JJA MAM SON

  ; -- time period of trends
  start_trend = 2021; 2021
  ende_trend  = 2050; 2050
  time_trend  = fspan(start_trend,ende_trend,ende_trend-start_trend+1)

  ; -- which ensemble member to decompose
  e = 15

  print(""+seas+" -- "+start_trend+"-"+ende_trend)
  print("  ensmem = "+e)

  ; -- grid
  grid = "native" ; native g025

  ; -- units of trend
  units = "\circC"


;--------------------------------------------------------------------------------

  function seas_and_trend(infield,seas,start,ende,start_trend,ende_trend)
  begin
    if seas .eq. "annual"
      tmp1 = month_to_annual(infield,1)
    else
      tmp1 = month_to_season(infield,seas)
    end if
    if units .eq. "%"
      tmp1 = ( tmp1 / conform_dims(dimsizes(tmp1),dim_avg_n_Wrap(tmp1(start_trend-start:ende_trend-start,:,:),0),(/1,2/)) )*100
    end if
    outfield = tofloat( regCoef_n(time_trend,tmp1(start_trend-start:ende_trend-start,:,:),0,0) * (ende_trend-start_trend) ) ; calculate trend
    delete(tmp1)
    copy_VarCoords(infield(0,:,:),outfield)
    return(outfield)
  end


  ; -- CESM LE data --
  start   = 1920
  ende    = 2100
  ensmem  = 40
  pathin = "/project/cas/flehner/LE/hist2rcp85/NA/"
  if recalc1 .eq. "yes"
    print("load dummy...")
    a = addfile(pathin+"tas_Amon_CESM1_hist2rcp85_r2i1p1_"+start+"01-"+ende+"12.nc","r")
    tmp     = a->tas
    dummy   = tofloat(tmp)
    copy_VarCoords(tmp,dummy)
    delete(tmp)
    dims    = dimsizes(dummy)
    printVarSummary(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)

    ; -- load data and create fields for 3x3 plot (a la Deser et al 2016, Fig. 9)
    tas_data    = new((/9,nlat,nlon/),float)
    a = addfile(pathin+"tas_Amon_CESM1_hist2rcp85_r"+e+"i1p1_"+start+"01-"+ende+"12.nc","r")
    tmp0  = a->tas
    tas_data(0,:,:) = seas_and_trend(tmp0,seas,start,ende,start_trend,ende_trend)
    a = addfile(pathin+"tas_Amon_CESM1_hist2rcp85_ensmean40_"+start+"01-"+ende+"12.nc","r")
    em    = a->tas
    tas_data(2,:,:) = seas_and_trend(em,seas,start,ende,start_trend,ende_trend)
    tas_data(1,:,:) = tas_data(0,:,:)-tas_data(2,:,:)

    a = addfile(pathin+"dynadj/tas_dynadj_Amon_CESM1_hist2rcp85_r"+e+"i1p1_"+start+"01-"+ende+"12.nc","r")
    tmp0    = a->tas
    tas_data(6,:,:) = seas_and_trend(tmp0,seas,start,ende,start_trend,ende_trend)
    a = addfile(pathin+"dynadj/tas_dynadj_Amon_CESM1_hist2rcp85_ensmean40_"+start+"01-"+ende+"12.nc","r")
    em_dynadj    = a->tas
    tas_data(8,:,:) = seas_and_trend(em_dynadj,seas,start,ende,start_trend,ende_trend)
    tas_data(7,:,:) = tas_data(6,:,:)-tas_data(8,:,:)

    tas_data(3,:,:) = tas_data(0,:,:)-tas_data(6,:,:)
    tas_data(4,:,:) = tas_data(1,:,:)-tas_data(7,:,:)
    tas_data(5,:,:) = tas_data(2,:,:)-tas_data(8,:,:)
    do i = 0,8
      copy_VarCoords(em(0,:,:),tas_data(i,:,:))
    end do

    psl_data    = new((/9,nlat,nlon/),float)
    a = addfile(pathin+"psl_Amon_CESM1_hist2rcp85_r"+e+"i1p1_"+start+"01-"+ende+"12.nc","r")
    tmp99 = a->PSL
    tmp0  = tmp99/100
    psl_data(0,:,:) = seas_and_trend(tmp0,seas,start,ende,start_trend,ende_trend)
    a = addfile(pathin+"psl_Amon_CESM1_hist2rcp85_ensmean30_"+start+"01-"+ende+"12.nc","r")
    tmp0  = a->PSL
    em    = tmp0/100
    delete(tmp0)
    psl_data(2,:,:) = seas_and_trend(em,seas,start,ende,start_trend,ende_trend)
    psl_data(1,:,:) = psl_data(0,:,:)-psl_data(2,:,:)

    a = addfile(pathin+"dynadj/psl_dynadj_Amon_CESM1_hist2rcp85_r"+e+"i1p1_"+start+"01-"+ende+"12.nc","r")
    tmp0    = rm_single_dims(a->psl)
    psl_data(3,:,:) = seas_and_trend(tmp0,seas,start,ende,start_trend,ende_trend)
    a = addfile(pathin+"dynadj/psl_dynadj_Amon_CESM1_hist2rcp85_ensmean30_"+start+"01-"+ende+"12.nc","r")
    em_dynadj    = a->psl
    psl_data(5,:,:) = seas_and_trend(em_dynadj,seas,start,ende,start_trend,ende_trend)
    psl_data(4,:,:) = psl_data(3,:,:)-psl_data(5,:,:)

    psl_data(6,:,:) = psl_data(0,:,:)-psl_data(3,:,:)
    psl_data(7,:,:) = psl_data(1,:,:)-psl_data(4,:,:)
    psl_data(8,:,:) = psl_data(2,:,:)-psl_data(5,:,:)
    do i = 0,8
      copy_VarCoords(tmp99(0,:,:),psl_data(i,:,:))
    end do
  else
    print_clock("load data...")
    ofile1 = pathin+"TREFHT_trends_"+seas+"_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"r")
    trends_model       = o->trends
  end if

;-----------------------------------------------------------------------------------------------









;====================================================================
; -- PLOTTING
print_clock("plotting...")

  wks_type = "eps" ; eps png
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500
     wks_type@wkHeight = 1500
  end if
  ofile = "/home/flehner/publication/postage_stamp_maps_random/"+get_script_prefix_name()+"_"+seas
  wks = gsn_open_wks(wks_type,ofile)

  gsn_define_colormap(wks,"amwg_blueyellowred")
  cmap = gsn_retrieve_colormap(wks)
  hsv = rgbhsv(cmap)
  hsv(10,1) = .47
  hsv(8,1) = .45    ; 1c+1d
  hsv(9,1) = .30    ; 1d
  hsv(4,:) = (/ (hsv(3,:)+hsv(5,:))/2.0 /)
  cmap2 = hsvrgb(hsv)
  gsn_define_colormap(wks,cmap2)
  newcolor = NhlNewColor(wks,0.42,0.42,0.42)
  nc0 = NhlNewColor(wks,0.95,0.95,0.95)
  nc1 = NhlNewColor(wks,0.7,0.7,0.7)
  nc3 = NhlNewColor(wks,0.5,0.5,0.5)
  nc4 = NhlNewColor(wks,0.30,0.30,0.30)
  nc2 = NhlNewColor(wks,139/256.,35/256.,35/256.)


  res = True
  res@mpFillOn = True
  res@mpFillColors = (/-1,0,-1,-1/)
  res@mpFillDrawOrder = "Draw"
  res@mpOutlineDrawOrder = "PostDraw"
  res@mpPerimDrawOrder = "PostDraw"
  res@mpOutlineDrawOrder = "Draw"
  res@cnLineDrawOrder = "PostDraw"
  res@cnFillDrawOrder = "PreDraw"
  res@mpGeophysicalLineColor = "gray50"
  res@mpUSStateLineColor = res@mpGeophysicalLineColor
  res@mpNationalLineColor = res@mpGeophysicalLineColor
  res@mpOutlineOn = True
  res@gsnAddCyclic = False
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpOutlineBoundarySets = "National" ;"GeophysicalAndUSStates" ; National USStates
  res@mpNationalLineColor = res@mpGeophysicalLineColor
  res@mpNationalLineThicknessF = res@mpGeophysicalLineThicknessF
     res@mpPerimOn    = True
     res@mpProjection = "Robinson"
     res@mpGridLatSpacingF =  90            ; change latitude  line spacing
     res@mpGridLonSpacingF = 180.           ; change longitude line spacing
     res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
     res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
     res@mpLimitMode = "LatLon"
     res@mpMinLonF = 202.
     res@mpMaxLonF = 303.
     res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2.0
     res@mpMinLatF = 21.5
     res@mpMaxLatF = 71.5

  res@gsnFrame     =  False
  res@gsnDraw = False
  res@cnLineLabelsOn = False
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@tiMainOn = False
  res@lbLabelBarOn = False
  res@gsnLeftStringOrthogonalPosF = -0.97
  res@gsnRightStringOrthogonalPosF = -0.97
  res@gsnLeftStringParallelPosF = .025
  res@gsnRightStringParallelPosF = .975
  res@txFontHeightF = 0.03

  res@cnLevelSelectionMode = "ExplicitLevels"
  ; -- for precip
  ; gsn_define_colormap(wks,"precip_diff_12lev")
  ; res@cnFillColors = (/2,3,4,5,7,9,10,11,12,13,14,15,16/)
  ; res@cnLevels = (/-40,-30,-20,-10,0,10,20,30,40/)
  ; -- for tas
  res@cnFillColors = (/3,4,5,7,9,10,11,12,13,14,15,16,17/)
  res@cnLevels = (/-4,-3,-2,-1,0,1,2,3,4/)

  res@gsnCenterString = ""
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tiMainOn = False

  res_psl = res
  delete(res_psl@cnLevels)
  res_psl@cnLevels        = fspan(-15,15,16)
  res_psl@gsnContourNegLineDashPattern = 1
  res_psl@cnFillOn        = False
  res_psl@cnLinesOn       = True
  res_psl@cnInfoLabelOn   = False


  sres = True
  sres@gsnDraw = False
  sres@gsnFrame = False
  sres@cnFillOn = False
  sres@cnLevelSelectionMode = "ExplicitLevels"
  dpat = "$_$_$_$_$_$_"
  n_i  = NhlNewDashPattern(wks,dpat)
  sres@gsnContourNegLineDashPattern = 14   ;14, 11, or n_i
  sres@cnLevels = fspan(-12.5,12.5,26)
  sres@cnLineThicknessF = 1.25
  sres@cnInfoLabelOn = False
  sres@tiMainOn = False
  sres@cnLineLabelsOn = False
  sres@gsnLeftString = ""
  sres@gsnRightString = ""
  sres@gsnCenterString = ""

  panres = True
  panres@gsnFrame = False
  panres@gsnMaximize = True
  panres@gsnPanelBottom = 0.05
  panres@gsnPanelTop = 0.97
  panres@txString = ""
  panres@txFontHeightF = 0.015
  panres@gsnPanelLabelBar = True
  panres@lbBoxLineColor = "gray20"
  panres@gsnPanelYWhiteSpacePercent = 2.
  panres@txFontHeightF = 0.020
  panres@pmLabelBarOrthogonalPosF = .01
  panres@lbLabelStride = 1
  panres@lbLabelFontHeightF = 0.012
  panres@pmLabelBarHeightF = 0.06
  panres@pmLabelBarWidthF = 0.55

  panres@lbTitlePosition          = "Bottom"
  panres@pmLabelBarOrthogonalPosF = -0.015 ; moves the label bar explicitly
  panres@lbLabelFontHeightF       = 0.012
  panres@lbTitleFontHeightF       = 0.012
  panres@lbBoxEndCapStyle = "TriangleBothEnds"
  panres@lbTitleString = seas+" temperature trend "+start_trend+"-"+ende_trend+" (~F34~0~F~C / "+dimsizes(time_trend)+" yr)"

  labels  = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

  ; ---------

  plots = new(9,graphic)
  plots_ov = plots
  res@mpPerimLineColor = "gray60"

  do i = 0,8
    res@gsnRightStringFontHeightF = 0.02
    res@gsnRightString = labels(i)
    plots(i) = gsn_csm_contour_map(wks,tas_data(i,:,:),res)
  end do
  do i = 0,8
    plots_ov(i) = gsn_csm_contour(wks,psl_data(i,:,:),res_psl)
    overlay(plots(i),plots_ov(i))
  end do

  gsres = True
  gsres@gsFillColor = "white"
  gsres@tfPolyDrawOrder = "PostDraw"

  dum = new(dimsizes(plots),"graphic")

  gsn_panel(wks,plots,(/3,3/),panres)

  frame(wks)
  delete(wks)

  system("convert -density 250x250 -crop 0x0 -trim +repage "+ofile+".eps "+ofile+".jpg")

  print_clock("...done")

end
