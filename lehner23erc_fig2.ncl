load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "~asphilli/ncl/asphilli_util.ncl"

; ---------------------------
; lehner23erc_fig2.ncl
; ---------------------------
; Code to produce Fig. 2 in Lehner and Deser (2023)
;
; Citation:
; Lehner, F., C. Deser (2023):
; Origin, importance, and predictive limits of internal climate variability
; Environmental Research: Climate, DOI: https://doi.org/10.1088/2752-5295/accf30
;
; Notes on code:
; - Only does calculations and plotting, not pre-processing
; - Requires concatenated monthly mean data from:
;   - CESM1 LE (b.e11.B20TRC5CNBDRD.f09_g16.*.cam.h0.*.nc and b.e11.BRCP85C5CNBDRD.f09_g16.*.cam.h0.*.nc)
; ---------------------------

begin

; -- parameters --
recalc1 = "yes" ; recalc CESM trends?
recalc2 = "yes" ; recalc obs trends?

; -- CESM LE or ME?
ens   = "le"; le me

; -- season
seas  = "DJF"; wateryear annual DJF JJA MAM SON

; -- time period of trends
start_trend = 2021; 2021
ende_trend  = 2050; 2050
time_trend  = fspan(start_trend,ende_trend,ende_trend-start_trend+1)

print(""+seas+" -- "+start_trend+"-"+ende_trend)

; -- grid
grid = "native" ; native g025

; -- units of trend
units = "\circC"


;--------------------------------------------------------------------------------
; -- CESM LE data --
if ens .eq. "le"
  ensmem  = 30;40
  start   = 1920
  ende    = 2100
  pathin = "/project/cas/flehner/LE/atm/TREFHT/"
else
  ensmem  = 15
  start   = 1920
  ende    = 2080
  pathin  = "/project/cas/flehner/ME/atm/TREFHT/"
end if
if recalc1 .eq. "yes"
  print("load dummy...")
  if grid .eq. "g025"
    a = addfile(pathin+"TREFHT.r2i1p1."+start+"01-"+ende+"12.g025.nc","r")
  else
    a = addfile(pathin+"TREFHT.r2i1p1."+start+"01-"+ende+"12.nc","r")
  end if
  tmp     = a->TREFHT
  dummy   = tofloat(tmp)
  copy_VarCoords(tmp,dummy)
  delete(tmp)
  dims    = dimsizes(dummy)
  printVarSummary(dummy)
  ntime   = dims(0)
  nlat    = dims(1)
  nlon    = dims(2)
  trends_model  = new((/ensmem,nlat,nlon/),float)
  do e = 1,ensmem
    print("  ensmem = "+e)
    if grid .eq. "g025"
      a = addfile(pathin+"TREFHT.r"+e+"i1p1."+start+"01-"+ende+"12.g025.nc","r")
    else
      a = addfile(pathin+"TREFHT.r"+e+"i1p1."+start+"01-"+ende+"12.nc","r")
    end if
    tmp0 = a->TREFHT
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    if units .eq. "%"
      tmp1 = ( tmp1 / conform_dims(dimsizes(tmp1),dim_avg_n_Wrap(tmp1(start_trend-start:ende_trend-start,:,:),0),(/1,2/)) )*100
    end if
    delete(tmp0)
    trends_model(e-1,:,:) = tofloat( regCoef_n(time_trend,tmp1(start_trend-start:ende_trend-start,:,:),0,0) * (ende_trend-start_trend) ) ; calculate trend
  end do
  delete(tmp1)
  copy_VarCoords(dummy(0,:,:),trends_model(0,:,:))
  print_clock("save data...")
  ofile1 = pathin+"TREFHT_trends_"+seas+"_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trends       = trends_model
else
  print_clock("load data...")
  ofile1 = pathin+"TREFHT_trends_"+seas+"_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  trends_model       = o->trends
end if
if units .eq. "mm"
  tmp = trends_model*85400*365 ; convert to mm/year
  copy_VarCoords(trends_model,tmp)
  delete(trends_model)
  trends_model = tmp
  delete(tmp)
end if
trends_model_ensmean  = dim_avg_n_Wrap(trends_model,0)
trends_model_stddev   = dim_stddev_n_Wrap(trends_model,0)
signal2noise          = trends_model_ensmean/trends_model_stddev
copy_VarCoords(trends_model_ensmean,signal2noise)

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
systemfunc("rm -rf "+ofile)
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
res@mpPerimDrawOrder = "PostDraw"
res@mpOutlineDrawOrder = "Draw"
res@cnLineDrawOrder = "PostDraw"
res@cnFillDrawOrder = "PreDraw"

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
; res@cnFillColors = (/3,4,5,7,9,10,11,12,13,14,15,16/)
res@cnFillColors = (/3,4,5,7,9,10,11,12,13,14,15,16,17/)
res@cnLevels = (/-4,-3,-2,-1,0,1,2,3,4/)

res@gsnCenterString = ""
res@gsnLeftString = ""
res@gsnRightString = ""
res@tiMainOn = False


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
panres@gsnPanelLabelBar = False;True
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


; ---------
res@lbBoxEndCapStyle = "TriangleBothEnds"
res@lbLabelFontHeightF = 0.03
res@lbTitleFontHeightF = 0.03

plot2 = new(ensmem+3,graphic)
res@mpPerimLineColor = "gray60"
do i = 0,ensmem-1; 0,29
   res@gsnRightString = (i+1)
   if i .eq. ensmem-1
     res@lbLabelBarOn = True
   else
     res@lbLabelBarOn = False
   end if
   plot2(i) = gsn_csm_contour_map(wks,trends_model(i,:,:),res)
end do

res_std = res
delete(res_std@cnLevels)
delete(res_std@cnFillColors)
res_std@cnFillColors = (/10,11,12,13,14,15,16,17/)
res_std@cnLevels = (/1,2,3,4/)
; res_std@lbLabelBarOn = True
res_std@lbLabelFontHeightF = 0.03
res_std@lbTitleFontHeightF = 0.03
res_std@lbBoxEndCapStyle = "TriangleHighEnd"
res_std@lbTitlePosition          = "Bottom"
res_std@lbTitleString = ""; "Standard deviation of trends [~F34~0~F~C / "+dimsizes(time_trend)+" yr]"
res_std@mpPerimLineColor = "gray20"
res_std@gsnRightString = ""; "Signal (EM)"

res_std@lbLabelBarOn = False
plot2(ensmem) = gsn_csm_contour_map(wks,trends_model_ensmean,res_std)

res_std@lbLabelBarOn = True
res_std@gsnRightString = ""; "Noise (STD)"
plot2(ensmem+1) = gsn_csm_contour_map(wks,trends_model_stddev,res_std)

res_std@lbLabelBarOn = False
res_std@gsnRightString = ""; "Signal/Noise (EM/STD)"
plot2(ensmem+2) = gsn_csm_contour_map(wks,signal2noise,res_std)



res@cnFillMode = "RasterFill"
delete(res@cnFillMode)
res@cnInfoLabelOn = False
res@gsnRightString = ""

gsres = True
gsres@gsFillColor = "white"
gsres@tfPolyDrawOrder = "PostDraw"

dum = new(dimsizes(plot2),"graphic")

panres@gsnPanelRowSpec = True
gsn_panel(wks,plot2,(/6,6,6,6,6,3/),panres)

frame(wks)
delete(wks)

system("convert -density 250x250 -crop 0x0 -trim +repage "+ofile+".eps "+ofile+".jpg")

print_clock("...done")

end
