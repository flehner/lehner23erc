load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; ---------------------------
; lehner23erc_fig1.ncl
; ---------------------------
; Code to produce Fig. 1 in Lehner and Deser (2023)
;
; Citation:
; Lehner, F., C. Deser (2023):
; Origin, importance, and predictive limits of internal climate variability
; Environmental Research: Climate, DOI: https://doi.org/10.1088/2752-5295/accf30
;
; Notes on code:
; - Only does calculations and plotting, not pre-processing
; - Requires concatenated monthly mean data from:
;   - CESM1 coupled piControl (b.e11.B1850C5CN.f09_g16.005.cam.h0.*.040001-220012)
;   - CESM1 uncoupled piControl (f.e11.F1850C5CN.f09_f09.001.cam.h0.*.000101-260012)
; ---------------------------

begin

system("echo  ")

; -- recalculate variance or load pre-calculated data
recalc1 = "y"; coupled piControl
recalc2 = "y"; uncoupled piControl

; -- parameters --
pathin    = "/glade/scratch/flehner/"
pathout   = "/glade/scratch/flehner/"
vari      = "TREFHT" ; PSL PRECT TREFHT
seas      = "annual" ; annual DJF
wl_list   = (/1,5,10,30,50/); running mean window length in years

do w = 0,dimsizes(wl_list)-1 ; loop over window length
  wl = wl_list(w); window length in years

  print(""+vari)
  print(""+seas)
  print(""+wl+" yr")

  ; -- coupled control --
  ifile = "b.e11.B1850C5CN.f09_g16.005.cam.h0."+vari+".040001-220012"
  if recalc1 .eq. "y"
    print_clock("load coupled control...")
    f1 = pathin+ifile+".nc"
    a = addfile(f1,"r")
    tmp       = a->$vari$
    ; -- sasonal mean
    if seas .eq. "annual"
      var = month_to_annual(tmp,1)
    else
      var = month_to_season(tmp,seas)
    end if
    ; -- running mean
    if wl .gt. 1
      coupled = runave_n_Wrap(var,wl,0,0)
    else
      coupled = var
    end if
    delete(var)
    ; -- calc variance
    coupled_var = dim_variance_n(coupled,0)
    copy_VarCoords(tmp(0,:,:),coupled_var)
    delete(tmp)
    ; -- save data
    ofile1 = pathin+ifile+"_variance_"+seas+"_wl"+wl+"yr"
    system("rm -rf "+ofile1+".nc")
    o = addfile(ofile1+".nc","c")
    o->$vari$   = coupled_var
    print_clock("...done")
  else
    ; -- load data
    ofile1 = pathin+ifile+"_variance_"+seas+"_wl"+wl+"yr"
    i = addfile(ofile1+".nc","r")
    coupled_var = i->$vari$
  end if

  ; -- uncoupled control --
  ifile = "f.e11.F1850C5CN.f09_f09.001.cam.h0."+vari+".000101-260012"
  if recalc2 .eq. "y"
    print_clock("load uncoupled control...")
    f1 = pathin+ifile+".nc"
    a = addfile(f1,"r")
    tmp       = a->$vari$
    ; -- sasonal mean
    if seas .eq. "annual"
      var = month_to_annual(tmp,1)
    else
      var = month_to_season(tmp,seas)
    end if
    ; -- running mean
    if wl .gt. 1
      uncoupled = runave_n_Wrap(var,wl,0,0)
    else
      uncoupled = var
    end if
    delete(var)
    ; -- calc variance
    uncoupled_var = dim_variance_n(uncoupled,0)
    copy_VarCoords(tmp(0,:,:),uncoupled_var)
    delete(tmp)
    ; -- save data
    ofile1 = pathin+ifile+"_variance_"+seas+"_wl"+wl+"yr"
    system("rm -rf "+ofile1+".nc")
    o = addfile(ofile1+".nc","c")
    o->$vari$   = uncoupled_var
    print_clock("...done")
  else
    ; -- load data
    ofile1 = pathin+ifile+"_variance_"+seas+"_wl"+wl+"yr"
    i = addfile(ofile1+".nc","r")
    uncoupled_var = i->$vari$
  end if

  ; -- calculate ratio
  ratio = uncoupled_var/coupled_var
  signi = ftest(uncoupled_var,2600/wl,coupled_var,1800/wl,0)
  copy_VarCoords(coupled_var,ratio)
  copy_VarCoords(coupled_var,signi)



  ; ===========================================================================
  ; PLOTTING

  minlat = -90
  maxlat = 90
  minlon = 0;-180
  maxlon = 360;180


  varp 		  = "/glade/u/home/flehner/" + "variance_ratio_"+vari+"_"+seas+"_wl"+wl+"yr"
  wks_type 	= "eps"
  wks 		  = gsn_open_wks(wks_type,varp)
  nbplot    = 1
  plot_t	  = new(nbplot,graphic)
  plot_p	  = new(nbplot,graphic)
  plot_ov	  = new(nbplot,graphic)

  ; Adam's colormap
  gsn_define_colormap(wks,"amwg_blueyellowred")
  cmap = gsn_retrieve_colormap(wks)
  hsv = rgbhsv(cmap)
  hsv(10,1) = .47
  hsv(8,1) = .45
  hsv(9,1) = .30
  hsv(4,:) = (/ (hsv(3,:)+hsv(5,:))/2.0 /)
  cmap2 = hsvrgb(hsv)
  gsn_define_colormap(wks,cmap2)
  newcolor = NhlNewColor(wks,0.42,0.42,0.42)
  nc0 = NhlNewColor(wks,0.95,0.95,0.95)
  nc1 = NhlNewColor(wks,0.7,0.7,0.7)
  nc3 = NhlNewColor(wks,0.5,0.5,0.5)
  nc4 = NhlNewColor(wks,0.30,0.30,0.30)
  nc2 = NhlNewColor(wks,139/256.,35/256.,35/256.)
  gsn_merge_colormaps(wks,"precip_diff_12lev",cmap2)


  res                     = True
  res@mpFillOn = True                       ; default = True
  res@mpFillColors = (/-1,-1,-1,0/)      ; map color fill ocean
  res@cnFillDrawOrder = "PreDraw"        ; color fill first
  res@mpOutlineDrawOrder = "Draw"      ; draw map outline after cnFill
  res@mpFillDrawOrder = "Draw"              ; draw map fill next
  res@mpPerimDrawOrder = "PostDraw"   ; draw map perimeter last
  res@gsnDraw              = False    ; do not draw picture
  res@gsnFrame             = False    ; do not advance frame
  res@gsnPaperOrientation  = "landscape"
  res@gsnAddCyclic = True
  res@mpLimitMode   = "LatLon"
  res@mpOutlineOn = True
  res@mpGeophysicalLineColor = "black";"gray42"
  res@mpCountyLineColor = "transparent"
  res@mpUSStateLineColor = "gray42"
  res@mpNationalLineColor = "black";"gray22"
  res@mpGeophysicalLineThicknessF = 1.4
  res@mpProjection            = "Robinson"

  res@mpPerimOn    = False
  res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter

  res@mpMinLatF            = minlat
  res@mpMaxLatF            = maxlat
  res@mpMinLonF            = minlon
  res@mpMaxLonF            = maxlon
  res@mpCenterLonF         = 10.;150.         ; defailt is 0 [GM]

  ; Define colors for filled contours
  res@cnFillOn             = True
  res@gsnSpreadColors      = True                  ; use full colormap
  res@cnLevelSelectionMode = "ExplicitLevels"

  res@cnLinesOn           = False         ; Turn lines off
  res@cnLineLabelsOn      = False         ; Turn labels off
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@mpGridLineColor      =  -1       ; Turn off grid lat/lon lines
  res@mpGridAndLimbOn      =  True
  res@lbLabelBarOn        = False
  res@cnMissingValFillPattern = "SolidFill"
  res@cnMissingValFillColor =-1
  res@tmXBOn               = False
  res@tmXTOn               = False
  res@tmYLOn               = False
  res@tmYROn               = False
  res@tiMainFontHeightF     = 0.015
  res@tiYAxisFontHeightF    = 0.015
  res@cnInfoLabelOn       = False

  res@cnLevels            = fspan(0.1,.9,9)
  res@cnFillColors        = (/18,19,20,21,22,23,24,25,26,27,28/)


  resov1       = True
  resov1@mpFillColors = (/-1,0,-1,0/)      ; map color fill ocean
  resov1@cnFillDrawOrder = "PreDraw"        ; color fill first
  resov1@mpOutlineDrawOrder = "Draw"      ; draw map outline after cnFill
  resov1@mpFillDrawOrder = "Draw"              ; draw map fill next
  resov1@mpPerimDrawOrder = "PostDraw"   ; draw map perimeter last
  resov1@gsnDraw              = False        ; don't draw yet
  resov1@gsnFrame             = False
  resov1@cnFillOn            = True
  resov1@cnLinesOn           = False
  resov1@cnLineLabelsOn      = False
  resov1@gsnLeftString = ""
  resov1@cnLineThicknessF      = 3
  resov1@cnMonoFillColor       = True
  resov1@cnMonoFillPattern     = False
  resov1@cnMonoFillScale       = True
  resov1@cnFillColor = "gray30"
  resov1@cnFillDotSizeF        = 0.003
  resov1@cnLevelSelectionMode  = "ExplicitLevels"
  resov1@cnFillPatterns        = (/-1,8/); 8=hatch, 6= crosshatch, 17=stippling
  resov1@cnFillScaleF          = 0.75
  resov1@cnLevels              = (/.05/)
  resov1@cnInfoLabelOn         = False
  resov1@lbLabelBarOn          = False


  ; ---------------------------------------------------------------------

  plot 	 = gsn_csm_contour_map(wks,ratio,res)
  plot_ov =  gsn_csm_contour(wks,signi,resov1)
  overlay(plot,plot_ov)

  ; ---------------------------------------------------------------------

  resP1                  = True                    ; modify the panel plot
  resP1@gsnMaximize     = True
  resP1@gsnPanelRowSpec = False;True
  resP1@gsnPanelBottom    = 0.1
  resP1@gsnPanelTop    = 1.0
  resP1@gsnPanelFigureStringsFontHeightF = 0.008

  resP1@txString                = ""
  resP1@lbTitlePosition         = "Bottom"
  resP1@lbTitleString           = "Fraction of variance explained by uncoupled CAM5 vs CESM1"

  resP1@gsnPanelYWhiteSpacePercent = 1
  resP1@gsnPanelXWhiteSpacePercent = 1
  resP1@pmLabelBarOrthogonalPosF = -.01
  resP1@lbBoxSizing       = "ExplicitSizing"
  resP1@lbOrientation     = "horizontal"



  resP1@gsnPanelFigureStringsPerimOn = False
  resP1@amJust          = "BottomLeft"

  resP1@gsnPaperOrientation = "portrait"
  resP1@txFontHeightF = 0.015
  resP1@pmLabelBarWidthF = 0.45
  resP1@pmLabelBarHeightF = 0.06
  resP1@lbLabelStride = 1
  resP1@lbBoxLineColor = "gray20"
  resP1@lbLabelFontHeightF  = 0.015
  resP1@lbTitleFontHeightF      = 0.015
  resP1@gsnPanelLabelBar    = True       ; common label bar
  resP1@gsnFrame            = False      ; don't advance frame yet
  resP1@lbLabelAutoStride   = True       ; nice label bar labels

  gsn_panel(wks,plot,(/1,1/),resP1)
  delete(wks)
  system("convert "+varp+".eps "+varp+".png")

end do ; end of wl loop


end
