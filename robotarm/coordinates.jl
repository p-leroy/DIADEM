export spherical2cart,cart2spherical,cart2elevation,elevation2cart,cart2elevation2,elevation2cart2

#coordinate manipulation functions

#-------------------------------------------------------------
# spherical to cartesian
#-------------------------------------------------------------
function spherical2cart(theta,phi)

x = cos.(deg2rad.(phi)).*sin.(deg2rad.(theta));
y = sin.(deg2rad.(phi)).*sin.(deg2rad.(theta));
z = (1 .+ 0 * cos.(deg2rad.(phi))) .* cos.(deg2rad.(theta));

return x,y,z
end
#-------------------------------------------------------------


#-------------------------------------------------------------
#cartesian to spherical
#-------------------------------------------------------------

function cart2spherical(x,y,z)

phi = rad2deg.(atan2.(y,x));
theta=rad2deg.(atan2.(x.*cos.(deg2rad.(phi))+y.*sin.(deg2rad.(phi)),z));

return theta, phi
end
#-------------------------------------------------------------

#-------------------------------------------------------------
#cartesian to elevation/azimuth
#-------------------------------------------------------------

function cart2elevation(x,y,z)

ele = rad2deg.(atan2.(y,z));
az = rad2deg.(atan2.(x,z.*cos.(deg2rad.(ele))+y.*sin.(deg2rad.(ele))));

return az,ele
end
#-------------------------------------------------------------

#-------------------------------------------------------------
# elevation/azimuth to cartesian
#-------------------------------------------------------------

function elevation2cart(az,ele)

x = sin.(deg2rad.(az));
y = sin.(deg2rad.(ele)).*cos.(deg2rad.(az));
z = cos.(deg2rad.(ele)).*cos.(deg2rad.(az));

return x,y,z
end
#-------------------------------------------------------------

#-------------------------------------------------------------
#cartesian to azimuth/elevation
#-------------------------------------------------------------

function cart2elevation2(x,y,z)

az = rad2deg.(atan2.(x,z));
ele = rad2deg.(atan2.(y,z.*cos.(deg2rad.(az))+y.*sin.(deg2rad.(az))));

return az,ele
end
#-------------------------------------------------------------


#-------------------------------------------------------------
# azimuth/elevation to cartesian
#-------------------------------------------------------------

function elevation2cart2(az,ele)

y = sin.(deg2rad.(ele));
x = sin.(deg2rad.(az)).*cos.(deg2rad.(ele));
z = cos.(deg2rad.(az)).*cos.(deg2rad.(ele));

return x,y,z
end
#---
