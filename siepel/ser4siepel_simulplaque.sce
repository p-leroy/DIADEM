deg2rad = %pi/180.0;

function [A]=steering_vector(k0,theta,phi,d,Matxy)
    xsource = d*sin(deg2rad*theta)*cos(deg2rad*phi);
    ysource = d*sin(deg2rad*theta)*sin(deg2rad*phi);
    zsource = d*cos(deg2rad*theta);
    Mat_dxy = abs(Matxy-(xsource+%i*ysource));
    Mat_d = abs(Mat_dxy+%i*zsource);
    Mat_theta = atan(Mat_dxy,zsource+0.0*Mat_dxy);
    A = exp(-%i*k0*Mat_d)./Mat_d;

endfunction


function [A]=steering_vector_planar(k0,x,y,z,Matxy)
    xsource = x;
    ysource = y;
    zsource = z;
    Mat_dxy = abs(Matxy-(xsource+%i*ysource));
    Mat_d = abs(Mat_dxy+%i*zsource);
    Mat_theta = atan(Mat_dxy,zsource+0.0*Mat_dxy);
    A = exp(-%i*k0*Mat_d)./Mat_d;

endfunction




// source location
d = 480000;
phi = 0.0;
theta = -10.0;
offset_source_x = 0.00;
xsource = d*sin(deg2rad*theta)*cos(deg2rad*phi)+offset_source_x;
ysource = d*sin(deg2rad*theta)*sin(deg2rad*phi);
zsource = d*cos(deg2rad*theta);





// source aperture : cos(theta).^n_source
n_source = 0;

// fréquence
f = 50;

// Dimensions de la plaque
Lx=600;
Ly=600;

dx = 0.5*300/f;
dy = 0.5*300/f;

nbx = ceil(Lx/dx)+1;
nby = ceil(Ly/dy)+1;

tabx = linspace(-0.5*Lx,0.5*Lx,nbx);
taby = linspace(-0.5*Ly,0.5*Ly,nby);

Matxy = (1.0+0.0*tabx')*tabx + %i*taby'*(1.0+0.0*taby);

k0 = 2*%pi*f/300.0;

Mat_dxy = abs(Matxy-(xsource+%i*ysource));
Mat_d = abs(Mat_dxy+%i*zsource);
Mat_theta = atan(Mat_dxy,zsource+0.0*Mat_dxy);

A_aperture = (cos(Mat_theta).^n_source).*steering_vector(k0,theta,phi,d,Matxy);



// observation
theta_obs = linspace(-30,30,3001);
phi_obs = 0.0;

E = zeros(length(theta_obs),1);

for i=1:length(theta_obs)
    A = steering_vector(k0,theta_obs(i),phi_obs,d,Matxy);
    E(i) = sum(A_aperture.*A);
end

//clf();
plot2d(theta_obs,20*log10(abs(E)));
xgrid();



