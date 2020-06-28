export compute_increment, compute_igloo, compute_igloo_v2, igloo_seq_file, igloo_seq_file_v2, read_sumup_igloo;

#######################################################################
## 			compute_increment
##
## Computation of the increment to be used for a balanced igloo ie
## an igloo with a phi increment enabling periodicity of the sampling
##
##	Input : angle (float)
##	Output : tab of acceptable increment
##
#######################################################################

function compute_increment(angle_max)

# factor(ContainerType, n::Integer) -> ContainerType
# Return the factorization of n stored in a ContainerType
tab_elem_init = [ 1; Primes.factor(Vector, round(Int, angle_max)) ];
tmp1 = [ 1 ];
tmp2 = collect( 0 : 1 : length(tab_elem_init) - 1);
for id = 0 : length(tab_elem_init) - 1
    tab_elem = tab_elem_init[ 1 .+ mod.(tmp2 .+ id, length(tab_elem_init)) ];
    for ie = 1 : length(tab_elem) - 1
        tmp = prod( tab_elem[1:ie] ) * tab_elem[ ie + 1 : length(tab_elem) ];
        tmp1 = [ tmp1; tmp ];
    end
end

tab_root = unique( tmp1[ findall( tmp1 .<= (0.5 * round(Int, angle_max) ) ) ] );
tab_div = [1 2 4 5 6 8 9 10 16 20 32 40 50 100];

tmp = [];
for ir = 1 : length(tab_root)
    cur_val = tab_root[ir] ./ tab_div;
    for jr = 1 : length(cur_val)
        if ( cur_val[ jr ] * 50 == round(Int, cur_val[ jr ] * 50) )
            tmp = [tmp; cur_val[jr]];
        end
    end
end

tab_increment = sort( unique( tmp ) );

return tab_increment;

end


#######################################################################
## 			compute_igloo
##
## Computation of the sampling for 0<=theta<=theta_max according to a
## given theta increment
##
##	Input : theta_step (float), theta_max (float)
##	Output : Matrix [Theta Phi]
##
#######################################################################

function compute_igloo(theta_step, theta_max)
tab_phi_step = compute_increment(360.0);

nb_elem = round(Int, theta_max / theta_step) + 1;
# OLD SYNTAX theta = collect(linspace(0.0, theta_max,nb_elem));
theta = collect(range(0.0, theta_max, length=nb_elem));

# Specific theta=0 location : phi_step = 5°
# OLD SYNTAX tmp = collect(linspace(0.0, 360.0, 73));
tmp = collect(range(0.0, 360.0, length=73));

igloo=[[] []];

for i = 1 : nb_elem

	phi_step = 0.0;
	if (abs(sin(deg2rad(theta[i]))) > 0.01)
		phi_step = abs(theta_step / sin(deg2rad((i-1) * theta_step)));
		phi_step = maximum(tab_phi_step[findall(tab_phi_step .<= phi_step)]);
	else
		phi_step = 60.0;
	end

	phi_shift = 0.5 * phi_step * isodd(i);
	nb_phi = round(Int, 360.0/phi_step) + 1;
	# OLD SYNTAX tab_phi = mod.(phi_shift+collect(Float64,linspace(0.0,360.0,nb_phi)),360.0);
    tab_phi = mod.(phi_shift .+ collect(Float64, range(0.0, 360.0, length=nb_phi)), 360.0);
    # OLD SYNTAX FOR round
    #igloo = [igloo; round.(theta_step*(i-1.0) .+ 0.0 * tab_phi[1:nb_phi-1], 3) round.(tab_phi[1:nb_phi-1], 3)];
	igloo = [igloo; round.(theta_step * (i-1.0) .+ 0.0 * tab_phi[1:nb_phi-1], digits=3) round.(tab_phi[1:nb_phi-1], digits=3)];

end

return igloo

@everywhere gc();

end

#######################################################################
## 			compute_igloo_v2
##
## Computation of the sampling for 0<=theta<=theta_max according to a
## given theta increment
##
##	Input : theta_step (float), theta_max (float)
##	Output : Matrix [Theta Phi]
##
#######################################################################

function compute_igloo_v2(theta_step, theta_max)

tab_increment_theta = compute_increment(theta_max);
tab_phi_step = compute_increment(360);

phi_step_min = 2 * theta_step; # on décale les phi une fois sur deux...

#phi_step_min = 2*rad2deg(asin(2*sin(0.5*deg2rad(theta_step))/sqrt(3)));

itmp = findmin(abs.(tab_increment_theta .- theta_step))[2];
theta_step = tab_increment_theta[itmp];
if (theta_step > phi_step_min) && (itmp > 1)
	theta_step = tab_increment_theta[itmp-1];
end

nb_elem = round(Int, theta_max/theta_step)+1;
theta = collect(range(0.0, theta_max, length=nb_elem));

#phi_step_min = 2*phi_step_min;
itmp = findmin(abs.(tab_phi_step .- phi_step_min))[2];
phi_step_min = tab_phi_step[itmp];
if (itmp > 1)
phi_step_min = tab_phi_step[itmp] > phi_step_min ? tab_phi_step[itmp-1] : tab_phi_step[itmp];
end

itmp = findall(mod.(tab_phi_step,phi_step_min) .== 0);
ntab_phi_step = tab_phi_step[itmp];

#itmp = indmax(tab_increment_theta[find(tab_increment_theta.<=theta_step)]);
#theta_step = tab_increment_theta[itmp];
#nb_elem = round(Int,theta_max/theta_step)+1;
#theta = collect(linspace(0.0,theta_max,nb_elem));

# Specific theta=0 location : phi_step = 5°
tmp = collect(range(0.0, 360.0, length=73));

igloo = [[] []];

for i = 1 : nb_elem

	phi_step = 0.0;
	if (abs(sin(deg2rad(theta[i]))) > 0.01)
		phi_step = abs(theta_step / sin(deg2rad((i-1) * theta_step)));
#		phi_step = 2*rad2deg(asin(2*sin(deg2rad(theta_step*0.5))/(sqrt(3)*sin(deg2rad(theta[i])))));
		itmp = findmin(abs.(ntab_phi_step .- phi_step))[2];
		phi_step = ntab_phi_step[itmp];
		if (i > 1)
		phi_step = ntab_phi_step[itmp]>phi_step ? ntab_phi_step[itmp-1] : ntab_phi_step[itmp];
		end

	else
		phi_step = 60.0;
	end

	phi_shift = 0.5 * phi_step * isodd(i);
	nb_phi = round(Int, 360.0 / phi_step) + 1;
	tab_phi = mod.(phi_shift .+ collect(Float64, range(0.0, 360.0, length=nb_phi)), 360.0);
	igloo = [igloo; round.(theta_step*(i-1.0) .+ 0.0 * tab_phi[1:nb_phi-1], digits=3) round.(tab_phi[1:nb_phi-1], digits=3)];

end

return igloo

@everywhere gc();

end

#######################################################################
## 			igloo_seq_file
##
## Generate sequential file used for measurement
##
##	Input : theta_max (float), azimuth_offset (float), roll_offset(float)
##		freq_min(float), freq_max(float), nb_freq (Int)
##		Antenna length(Float), Antenna heigth(Float), Antenna width(Float)
##		facility (string) : CAMILL or CACENDRA
##		path of the directory where to save the output files on the current computer,
##		path of the directory where to locate files on the controler computer,
##		Prefixe of the measured file
##	Output : Sequential file : SEQ.seq,
##		 (theta,phi) data sample location : NURBS.txt
##
#######################################################################

function igloo_seq_file(theta_max,freq_min,freq_max,nb_freq,
			AUT_depth,AUT_heigh,AUT_width,
			facility,
			PATH_FILE,
			PATH_SEQFILE;
			polar_angles=[0.0 90.0],
			azimuth_offset=0.0,
			roll_offset=0.0,
			DataFiles_prefixe="EvEh",
			ProjectName="3D",
			Operator="LLC",
			flag_half_sphere = 1,
			tdrift_periodicity = 450,
			freq_list=[])
cd();

if !(isdir(PATH_FILE))
	cd();
	mkdir(PATH_FILE);
	cd();
	cd(PATH_FILE);
else
	cd(PATH_FILE);
	if isfile("SEQ.seq") rm("SEQ.seq"); end;
	if isfile("sumup_igloo.txt") rm("sumup_igloo.txt"); end;
	if isfile("igloo.txt") rm("igloo.txt"); end;

	i = 0
	tmp = readdir(PATH_FILE)
	#tutu = [ismatch(r"config",i) for i in tmp]
    tutu = [occursin(r"config", i) for i in tmp]

	if (tutu!=[])
	itutu = findall(tutu);
	i = 0;

	for i in itutu
#		print(i)
		rm(tmp[i]);
	end
	end
end;

PATH_SEQFILE = join(split(PATH_SEQFILE, "/"), "\\")
DataFiles_dirgene = PATH_SEQFILE
ProjectName = "3D"
Operator = "LLC"

ScanAxis = 0
MOVE_Axis = 0
StepAxis = 0

if (facility.=="CAMILL_CATR")
	ScanAxis = 3 ; # roulis pour la chambre cm
	acc_scan = 8; # accélération de l'axe du scan
	v0_scan = 2;  # vitesse de l'axe du scan
	t0_scan = v0_scan / acc_scan;

	MOVE_Axis = 1 ; # axe azimuth
	acc_move = 8; # accélération de l'axe move
	v0_move = 2;  # vitesse de l'axe move
	t0_move = v0_move / acc_move;

	StepAxis= 6;
	acc_polar = 10; # accélération de l'axe polar
	v0_polar = 10;  # vitesse de l'axe polar

elseif (facility.=="CAMILL_FF")
	ScanAxis = 3 ; # roulis pour la chambre cm
	acc_scan = 8; # accélération de l'axe du scan
	v0_scan = 2;  # vitesse de l'axe du scan
	t0_scan = v0_scan / acc_scan;

	MOVE_Axis = 1 ; # axe azimuth
	acc_move = 8; # accélération de l'axe move
	v0_move = 2;  # vitesse de l'axe move
	t0_move = v0_move / acc_move;

	StepAxis= 5;
	acc_polar = 8; # accélération de l'axe polar
	v0_polar = 2;  # vitesse de l'axe polar

elseif (facility.=="CACENDRA")
	ScanAxis = 3 ; # roulis pour la chambre cm
	acc_scan = 8; # accélération de l'axe du scan
	v0_scan = 2;  # vitesse de l'axe du scan
	t0_scan = v0_scan/acc_scan;

	MOVE_Axis = 1 ; # axe azimuth
	acc_move = 8; # accélération de l'axe move
	v0_move = 2;  # vitesse de l'axe move
	t0_move = v0_move/acc_move;

	StepAxis = 4;
	acc_polar = 8; # accélération de l'axe polar
	v0_polar = 2;  # vitesse de l'axe polar

else
	println("Facility name is not valid or specified. Current facilities : CAMILL_CATR; CAMILL_FF ; CACENDRA")
	return -1
end

tab_increment_theta = compute_increment(theta_max);
tab_increment_phi = compute_increment(360);

# Theta step computation
Radius = sqrt((AUT_width*0.5)^2 + (AUT_heigh*0.5)^2 + (AUT_depth*0.5)^2);

Nmodes = 2*pi*Radius*freq_max/300.0;
delta_theta = 180.0/Nmodes;

if (flag_half_sphere!=0)
    delta_theta = delta_theta * sqrt(2);
    theta_max = min(90, theta_max);
end

#j = indmax(tab_increment_theta[findall(tab_increment_theta .<= delta_theta)]);
j = findmax(tab_increment_theta[findall(tab_increment_theta .<= delta_theta)])[2];

theta_step = tab_increment_theta[j];

flag_shift = 1;

# controle de la dérive en température
itemp = tdrift_periodicity;

# definition de la plage de fréquence en GHz

tab_freq = collect(range(freq_min, freq_max, length=nb_freq));

# computation of the samples locations
igloo = compute_igloo(theta_step, theta_max);
phi_max = maximum(igloo[:,2]);
phi_min = minimum(igloo[:,2]);
delta_phi_RJ = min(phi_min, 360-phi_max);
if (delta_phi_RJ == 0.0)
	delta_phi_RJ = max(phi_min, 360-phi_max);
end

igloo[:,2] = mod.(igloo[:,2] .+ roll_offset,360);

if (facility.=="CAMILL_CATR") || (facility.=="CAMILL_FF")

	tmp = findall(igloo[:,2] .> 180.0);
	igloo[tmp,1] = -1.0 * igloo[tmp,1];
	igloo[tmp,2] = -180.0 + igloo[tmp,2];

	# For CAMILL roll axis has to be negative...
#	igloo[tmp,:] = -1*igloo[tmp,:];
#	igloo[:,1]=abs.(sign.(igloo[:,1])).*igloo[:,1];
#	igloo[:,2]=abs.(sign.(igloo[:,2])).*igloo[:,2];

	igloo[:,1] = -1.0 * igloo[:,1];
	igloo[:,2] = -180.0 + igloo[:,2];
end

igloo[:,1] = igloo[:,1] .+ azimuth_offset;

tmp = ceil.(floor.(igloo[:,1], digits=4) .- 0.0001, digits=3) + im*ceil.(floor.(igloo[:,2], digits=4) .- 0.0001, digits=3);
tmp = unique(tmp);

igloo = [];
GC.gc();
igloo = zeros(Float64, size(tmp)[1], 2);
igloo[:,1] = real.(tmp);
igloo[:,2] = imag.(tmp);

#igloo = sortrows(igloo);
igloo = sortslices(igloo, dims=1);

tab_theta = unique(igloo[:,1]);
nb_theta = length(tab_theta);

# Sequential file head part

# OLD SYNTAX date = Date();
date = Dates.Date(1);
entete = @sprintf("Ant32 RasterScan Setup File V3.0\nDate=%s\nData#1\nProject=%s\nOperator=%s\nRemarks1=\nRemarks2= \nRemarks3= \nXyzImgAxis=0\nXyz ImgLabel= 0\nScanAxis=%d\n",date,ProjectName,Operator,ScanAxis);

# definition de la partie polarisation
ScanSpeed = 100;
ScanMode = "StepMode";
ScanMotionMode = "BiDirectional";
AxisLatency = 10;

Nb_polar = length(polar_angles);

Polar1 = polar_angles[1];
if Nb_polar > 1
	Polar2 = polar_angles[2];
end;

str_polar = [];
if (Nb_polar == 1)
    str_polar = @sprintf("ScanSpeed= %d\nScanMode= %s\nScanMotionMode= %s\nScanFreeRotation= False\nAxisLatency= %d\nStepAxis= %d\nMeasParam#=0\nStep#= %d\nStep 1= %f\n",ScanSpeed,ScanMode,ScanMotionMode,AxisLatency,StepAxis,Nb_polar,Polar1);
else
    str_polar = @sprintf("ScanSpeed= %d\nScanMode= %s\nScanMotionMode= %s\nScanFreeRotation= False\nAxisLatency= %d\nStepAxis= %d\nMeasParam#=0\nStep#= %d\nStep 1= %f\nStep 2= %f\n",ScanSpeed,ScanMode,ScanMotionMode,AxisLatency,StepAxis,Nb_polar,Polar1,Polar2);
end;

# definition de la plage de fréquence
str_freq = "";

if freq_list!=[]
# Liste de fréquence
nb_freq = length(freq_list);
str_freq = @sprintf("Freq#=%d\n", nb_freq);
for i=1:nb_freq,
  tmp = @sprintf("Freq %d=%8.6f\n", i, freq_list[i]);
  str_freq = string(str_freq,tmp);
end;

else
# Sweep frequentiel
str_freq = @sprintf("SweepStart=%8.6f\nSweepStop=%8.6f\nSweepPoints=%d\n", freq_min, freq_max, nb_freq);
end

# Sequence definition and file creation
index = 1;
nbdata = 0;

# First measurement : special case ie azimuth 0 check for Rotary joint characterization
if (index == 1)
	ptf = open("SEQ.seq","w");

	@printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=FWD\n",0,MOVE_Axis,-1);
	close(ptf);

	delta_phi = delta_phi_RJ;

	#j = indmin(abs.(tab_increment_phi .- delta_phi));
    j = findmin(abs.(tab_increment_phi .- delta_phi))[2];
	ScanInc = tab_increment_phi[j];

	if (facility.=="CAMILL_CATR")||(facility.=="CAMILL_FF")
		ScanStart = -180;
		ScanStop = 180 + ScanStart;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	else
		ScanStart = 0;
		ScanStop = 360 + ScanStart - ScanInc;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	end

	DataFiles=string(DataFiles_dirgene, DataFiles_prefixe, @sprintf("_scan%d.amp", index));
	RSSFile = @sprintf("config%d.rss", index);
	ptf = open(RSSFile, "w");
	@printf(ptf, "%s", entete);
	@printf(ptf, "ScanStart= %f\nScanStop= %f\nScanInc= %f\n", ScanStart, ScanStop, ScanInc);
	@printf(ptf, "%s", str_polar);
	@printf(ptf, "%s", str_freq);
	@printf(ptf, "DataFiles= %s\n#End", DataFiles);
	close(ptf);

	ptf = open("SEQ.seq", "a");
	@printf(ptf, "[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=REV\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n", 2*(index-1)+1, MOVE_Axis, azimuth_offset, 2*index, DataFiles_dirgene, RSSFile);
	close(ptf);

	index = index + 1;

end;

for i=1:nb_theta

  if (mod(i, itemp) == 0)
# Special case : reference location for thermal drift estimation
	delta_phi = delta_phi_RJ;

	j = indmin(abs(tab_increment_phi - delta_phi));
	ScanInc = tab_increment_phi[j];

	if (facility.=="CAMILL_CATR")||(facility.=="CAMILL_FF")
		ScanStart = -180;
		ScanStop = 180 + ScanStart;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	else
		ScanStart = 0;
		ScanStop = 360 + ScanStart - ScanInc;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	end

	DataFiles=string(DataFiles_dirgene,DataFiles_prefixe, @sprintf("_scan%d.amp", index));
	RSSFile = @sprintf("config%d.rss", index);
	ptf = open(RSSFile, "w");
	@printf(ptf, "%s",entete);
	@printf(ptf, "ScanStart= %f\nScanStop= %f\nScanInc= %f\n", ScanStart, ScanStop, ScanInc);
	@printf(ptf, "%s", str_polar);
	@printf(ptf, "%s", str_freq);
	@printf(ptf, "DataFiles= %s\n#End", DataFiles);
	close(ptf);

	ptf = open("SEQ.seq","a");
	@printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=REV\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n",2*(index-1)+1,MOVE_Axis,azimuth_offset,2*index,DataFiles_dirgene,RSSFile);
	close(ptf);

	index = index + 1;

  end

  theta_current = tab_theta[i];
  tmp = findall(abs.(igloo[:,1] .- theta_current).<=0.01);
  phi_current = igloo[tmp,2];
  ScanInc = phi_current[2]-phi_current[1];

  ScanStart = phi_current[1];
  ScanStop = phi_current[length(tmp)];

  nb_phi =  round(Int, (ScanStop-ScanStart)/ScanInc + 1);
  ScanInc = round((ScanStop-ScanStart)/(nb_phi-1), digits=3);

  nbdata = nbdata + nb_phi;

  DataFiles=string(DataFiles_dirgene,DataFiles_prefixe,@sprintf("_scan%d.amp",index));
  RSSFile = @sprintf("config%d.rss",index);
  ptf = open(RSSFile,"w");
  @printf(ptf, "%s",entete);
  @printf(ptf, "ScanStart= %f\nScanStop= %f\nScanInc= %f\n", ScanStart, ScanStop, ScanInc);
  @printf(ptf, "%s", str_polar);
  @printf(ptf, "%s", str_freq);
  @printf(ptf, "DataFiles= %s\n#End", DataFiles);
  close(ptf);

  ptf = open("SEQ.seq","a");
  @printf(ptf, "[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=FWD\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n", 2*(index-1)+1, MOVE_Axis, theta_current, 2*index, DataFiles_dirgene, RSSFile);
  close(ptf);
  index = index + 1;
end

 delta_phi = delta_phi_RJ;

	j = findmin(abs.(tab_increment_phi .- delta_phi))[2];
	ScanInc = tab_increment_phi[j];

	if (facility.=="CAMILL_CATR")||(facility.=="CAMILL_FF")
		ScanStart = -180;
		ScanStop = 180 + ScanStart;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	else
		ScanStart = 0;
		ScanStop = 360 + ScanStart - ScanInc;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	end

	DataFiles=string(DataFiles_dirgene,DataFiles_prefixe,@sprintf("_scan%d.amp",index));
	RSSFile = @sprintf("config%d.rss",index);
	ptf = open(RSSFile,"w");
	@printf(ptf, "%s", entete);
	@printf(ptf, "ScanStart= %f\nScanStop= %f\nScanInc= %f\n", ScanStart, ScanStop, ScanInc);
	@printf(ptf, "%s", str_polar);
	@printf(ptf, "%s", str_freq);
	@printf(ptf, "DataFiles= %s\n#End", DataFiles);
	close(ptf);

	ptf = open("SEQ.seq","a");
	@printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=REV\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n", 2*(index-1)+1, MOVE_Axis, azimuth_offset, 2*index, DataFiles_dirgene, RSSFile);
	close(ptf);

	index = index + 1;

writedlm("igloo.txt", igloo);

#/ Transmission measurement estimation
tsweep = 1.0;
tdata_transfer = 0.2;

#// ANT32 config file loading time
tconfig_file = 2;

accelaration_phi = acc_scan;
vitesse_phi = v0_scan;
acceleration_polar = acc_polar;
vitesse_polar = v0_polar;
acceleration_theta = acc_move;
vitesse_theta = v0_move;

tvnom_phi = vitesse_phi/accelaration_phi;

tmp = igloo[:,1];

theta= unique(tmp);
dtheta = abs(theta[2]-theta[1]);
t0 = vitesse_theta/acceleration_theta;
d0 = 0.5*acceleration_theta*t0^2;
tdtheta = 2*t0 + (dtheta-2*d0)/vitesse_theta;

dpolar = 90;
t0 = vitesse_polar/acceleration_polar;
d0 = 0.5*acceleration_polar*t0^2;
tdpolar = 2*t0 + (dpolar-2*d0)/vitesse_polar;

estimated_time = -tdtheta;
for itheta=1:length(theta)
    itmp = findall(abs.(igloo[:,1] .- theta[itheta]) .< 0.001);
    phi = igloo[itmp,2];

    dphi = abs(phi[2]-phi[1]);
    nbphi = length(phi)-1;

    t0 = vitesse_phi/accelaration_phi;
    d0 = 0.5*accelaration_phi*t0^2;
    tdphi = 2*t0 + (dphi-2*d0)/vitesse_phi;

    estimated_time = estimated_time + (nbphi*(tdphi+tsweep) + tdpolar + tdata_transfer)*2+tdtheta+tconfig_file;

end

s = mod(estimated_time,60);
m = mod((estimated_time-s)/60,60);
h = (estimated_time-s-60*m)/3600;

@printf("\nEstimated time = %dh %2.2dm %2.2ds\n",h,m,s);

end

#######################################################################
## 			igloo_seq_file_v2
##
## Generate sequential file used for measurement
##
##	Input : theta_max (float), azimuth_offset (float), roll_offset(float)
##		freq_min(float), freq_max(float), nb_freq (Int)
##		Antenna length(Float), Antenna heigth(Float), Antenna width(Float)
##		facility (string) : CAMILL or CACENDRA
##		path of the directory where to save the output files on the current computer,
##		path of the directory where to locate files on the controler computer,
##		Prefixe of the measured file
##	Output : Sequential file : SEQ.seq,
##		 (theta,phi) data sample location : NURBS.txt
##
#######################################################################


function igloo_seq_file_v2(theta_max,freq_min,freq_max,nb_freq,
			AUT_depth,AUT_heigh,AUT_width,
			facility,
			PATH_FILE,
			PATH_SEQFILE;
			polar_angles=[0.0 90.0],
			azimuth_offset=0.0,
			roll_offset=0.0,
			DataFiles_prefixe="EvEh",
			ProjectName="3D",
			Operator="LLC",
			flag_half_sphere = 1,
			tdrift_periodicity = 450,
			freq_list=[])

cd();

if !(isdir(PATH_FILE))
	cd();
	mkdir(PATH_FILE);
	cd();
	cd(PATH_FILE);
else
	cd(PATH_FILE);
	if isfile("SEQ.seq")	rm("SEQ.seq"); end;
	if isfile("sumup_igloo.txt") rm("sumup_igloo.txt"); end;
	if isfile("igloo.txt") rm("igloo.txt"); end;

	i=0;
	tmp = readdir(PATH_FILE);
	tutu=[ismatch(r"config",i) for i in tmp];
	if (tutu!=[])
	itutu = find(tutu);
	i = 0;

	for i in itutu
#		print(i)
		rm(tmp[i]);
	end
	end
end;

PATH_SEQFILE=join(split(PATH_SEQFILE,"/"),"\\");
DataFiles_dirgene = PATH_SEQFILE;
ProjectName = "3D";
Operator = "LLC";
ScanAxis = 0;
MOVE_Axis = 0;
StepAxis = 0;

if (facility.=="CAMILL_CATR")
	ScanAxis = 3 ; # roulis pour la chambre cm
	acc_scan = 8; # accélération de l'axe du scan
	v0_scan = 2;  # vitesse de l'axe du scan
	t0_scan = v0_scan/acc_scan;

	MOVE_Axis = 1 ; # axe azimuth
	acc_move = 8; # accélération de l'axe move
	v0_move = 2;  # vitesse de l'axe move
	t0_move = v0_move/acc_move;

	StepAxis= 6;
	acc_polar = 10; # accélération de l'axe polar
	v0_polar = 10;  # vitesse de l'axe polar

elseif (facility.=="CAMILL_FF")
	ScanAxis = 3 ; # roulis pour la chambre cm
	acc_scan = 8; # accélération de l'axe du scan
	v0_scan = 2;  # vitesse de l'axe du scan
	t0_scan = v0_scan/acc_scan;

	MOVE_Axis = 1 ; # axe azimuth
	acc_move = 8; # accélération de l'axe move
	v0_move = 2;  # vitesse de l'axe move
	t0_move = v0_move/acc_move;

	StepAxis= 5;
	acc_polar = 8; # accélération de l'axe polar
	v0_polar = 2;  # vitesse de l'axe polar

elseif (facility.=="CACENDRA")
	ScanAxis = 3 ; # roulis pour la chambre cm
	acc_scan = 8; # accélération de l'axe du scan
	v0_scan = 2;  # vitesse de l'axe du scan
	t0_scan = v0_scan/acc_scan;

	MOVE_Axis = 1 ; # axe azimuth
	acc_move = 8; # accélération de l'axe move
	v0_move = 2;  # vitesse de l'axe move
	t0_move = v0_move/acc_move;

	StepAxis= 4;
	acc_polar = 8; # accélération de l'axe polar
	v0_polar = 2;  # vitesse de l'axe polar

else
	println("Facility name is not valid or specified. Current facilities : CAMILL_CATR; CAMILL_FF ; CACENDRA")
	return -1
end

tab_increment_theta = compute_increment(theta_max);
tab_increment_phi = compute_increment(360);

# Theta step computation
Radius = sqrt((AUT_width*0.5)^2+(AUT_heigh*0.5)^2+(AUT_depth*0.5)^2);

Nmodes = 2*pi*Radius*freq_max/300.0;

delta_theta = 180.0/Nmodes;

if (flag_half_sphere!=0)
#    delta_theta = delta_theta*sqrt(2);
    theta_max = min(90,theta_max);
end


j=indmax(tab_increment_theta[find(tab_increment_theta.<=delta_theta)]);
theta_step=tab_increment_theta[j];

#valeur utilisée dans igloo2
delta_theta = rad2deg(2.0*asin(sin(deg2rad(delta_theta*0.5))*sqrt(3)*0.5));


flag_shift = 1;

# controle de la dérive en température
itemp = tdrift_periodicity;

# definition de la plage de fréquence en GHz

tab_freq=collect(linspace(freq_min,freq_max,nb_freq));

# computation of the samples locations
#igloo = compute_igloo_v2(delta_theta, theta_max);

igloo = compute_igloo_v2(theta_step, theta_max);

phi_max = maximum(igloo[:,2]);
phi_min = minimum(igloo[:,2]);
delta_phi_RJ = min(phi_min,360-phi_max)*0.5;
if (delta_phi_RJ==0.0)
	delta_phi_RJ = max(phi_min,360-phi_max)*0.5;
end

igloo[:,2] = mod.(igloo[:,2]+roll_offset,360);

if (facility.=="CAMILL_CATR")||(facility.=="CAMILL_FF")


	tmp = find(igloo[:,2].>180.0);
	igloo[tmp,1]=-1.0*igloo[tmp,1];
	igloo[tmp,2]=-180.0+igloo[tmp,2];

	# For CAMILL roll axis has to be negative...
#	igloo[tmp,:] = -1*igloo[tmp,:];
#	igloo[:,1]=abs.(sign.(igloo[:,1])).*igloo[:,1];
#	igloo[:,2]=abs.(sign.(igloo[:,2])).*igloo[:,2];

	igloo[:,1]=-1.0*igloo[:,1];
	igloo[:,2]=-180.0+igloo[:,2];
end

igloo[:,1] = igloo[:,1] + azimuth_offset;

tmp = ceil.(floor.(igloo[:,1],4)-0.0001,3)+im*ceil.(floor.(igloo[:,2],4)-0.0001,3);
tmp= unique(tmp);

igloo=[];
gc();
igloo=zeros(Float64,size(tmp)[1],2);
igloo[:,1]=real.(tmp);
igloo[:,2]=imag.(tmp);

igloo = sortrows(igloo);

tab_theta = unique(igloo[:,1]);
nb_theta = length(tab_theta);

# Sequential file head part

date = Date();

entete = @sprintf("Ant32 RasterScan Setup File V3.0\nDate=%s\nData#1\nProject=%s\nOperator=%s\nRemarks1=\nRemarks2= \nRemarks3= \nXyzImgAxis=0\nXyz ImgLabel= 0\nScanAxis=%d\n",date,ProjectName,Operator,ScanAxis);



# definition de la partie polarisation
ScanSpeed= 100;
ScanMode= "StepMode";
ScanMotionMode= "BiDirectional";
AxisLatency= 10;

Nb_polar = length(polar_angles);

Polar1 = polar_angles[1];
if Nb_polar>1
	Polar2 = polar_angles[2];
end;

str_polar=[];
if (Nb_polar==1)
    str_polar=@sprintf("ScanSpeed= %d\nScanMode= %s\nScanMotionMode= %s\nScanFreeRotation= False\nAxisLatency= %d\nStepAxis= %d\nMeasParam#=0\nStep#= %d\nStep 1= %f\n",ScanSpeed,ScanMode,ScanMotionMode,AxisLatency,StepAxis,Nb_polar,Polar1);
else
    str_polar=@sprintf("ScanSpeed= %d\nScanMode= %s\nScanMotionMode= %s\nScanFreeRotation= False\nAxisLatency= %d\nStepAxis= %d\nMeasParam#=0\nStep#= %d\nStep 1= %f\nStep 2= %f\n",ScanSpeed,ScanMode,ScanMotionMode,AxisLatency,StepAxis,Nb_polar,Polar1,Polar2);
end;






# definition de la plage de fréquence
str_freq="";

if freq_list!=[]
# Liste de fréquence
nb_freq = length(freq_list);
str_freq=@sprintf("Freq#=%d\n",nb_freq);
for i=1:nb_freq,
  tmp = @sprintf("Freq %d=%8.6f\n",i,freq_list[i]);
  str_freq=string(str_freq,tmp);
end;

else
# Sweep frequentiel
str_freq=@sprintf("SweepStart=%8.6f\nSweepStop=%8.6f\nSweepPoints=%d\n",freq_min,freq_max,nb_freq);
end


# Sequence definition and file creation
index = 1;
nbdata=0;


# First measurement : special case ie azimuth 0 check for Rotary joint characterization
if (index==1)
	ptf = open("SEQ.seq","w");

	@printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=FWD\n",0,MOVE_Axis,-1);
	close(ptf);


	delta_phi = delta_phi_RJ;

	j=indmin(abs.(tab_increment_phi-delta_phi));
	ScanInc = tab_increment_phi[j];

	if (facility.=="CAMILL_CATR")||(facility.=="CAMILL_FF")
		ScanStart = -180;
		ScanStop = 180 + ScanStart;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	else
		ScanStart = 0;
		ScanStop = 360 + ScanStart - ScanInc;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	end

	DataFiles=string(DataFiles_dirgene,DataFiles_prefixe,@sprintf("_scan%d.amp",index));
	RSSFile = @sprintf("config%d.rss",index);
	ptf = open(RSSFile,"w");
	@printf(ptf,"%s",entete);
	@printf(ptf,"ScanStart= %f\nScanStop= %f\nScanInc= %f\n",ScanStart, ScanStop, ScanInc);
	@printf(ptf,"%s",str_polar);
	@printf(ptf,"%s",str_freq);
	@printf(ptf,"DataFiles= %s\n#End",DataFiles);
	close(ptf);

	ptf = open("SEQ.seq","a");
	@printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=REV\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n",2*(index-1)+1,MOVE_Axis,azimuth_offset,2*index,DataFiles_dirgene,RSSFile);
	close(ptf);

	index = index + 1;

end;


for i=1:nb_theta

  if (mod(i,itemp)==0)
# Special case : reference location for thermal drift estimation
	delta_phi = delta_phi_RJ;

	j=indmin(abs(tab_increment_phi-delta_phi));
	ScanInc = tab_increment_phi[j];

	if (facility.=="CAMILL_CATR")||(facility.=="CAMILL_FF")
		ScanStart = -180;
		ScanStop = 180 + ScanStart;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	else
		ScanStart = 0;
		ScanStop = 360 + ScanStart - ScanInc;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	end

	DataFiles=string(DataFiles_dirgene,DataFiles_prefixe,@sprintf("_scan%d.amp",index));
	RSSFile = @sprintf("config%d.rss",index);
	ptf = open(RSSFile,"w");
	@printf(ptf,"%s",entete);
	@printf(ptf,"ScanStart= %f\nScanStop= %f\nScanInc= %f\n",ScanStart, ScanStop, ScanInc);
	@printf(ptf,"%s",str_polar);
	@printf(ptf,"%s",str_freq);
	@printf(ptf,"DataFiles= %s\n#End",DataFiles);
	close(ptf);

	ptf = open("SEQ.seq","a");
	@printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=REV\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n",2*(index-1)+1,MOVE_Axis,azimuth_offset,2*index,DataFiles_dirgene,RSSFile);
	close(ptf);

	index = index + 1;

  end



  theta_current = tab_theta[i];
  tmp = find(abs.(igloo[:,1]-theta_current).<=0.01);
  phi_current = igloo[tmp,2];
  ScanInc = phi_current[2]-phi_current[1];

  ScanStart = phi_current[1];
  ScanStop = phi_current[length(tmp)];

  nb_phi =  round(Int,(ScanStop-ScanStart)/ScanInc + 1);
  ScanInc = round((ScanStop-ScanStart)/(nb_phi-1),3);


  nbdata = nbdata + nb_phi;

  DataFiles=string(DataFiles_dirgene,DataFiles_prefixe,@sprintf("_scan%d.amp",index));
  RSSFile = @sprintf("config%d.rss",index);
  ptf = open(RSSFile,"w");
  @printf(ptf,"%s",entete);
  @printf(ptf,"ScanStart= %f\nScanStop= %f\nScanInc= %f\n",ScanStart, ScanStop, ScanInc);
  @printf(ptf,"%s",str_polar);
  @printf(ptf,"%s",str_freq);
  @printf(ptf,"DataFiles= %s\n#End",DataFiles);
  close(ptf);

  ptf = open("SEQ.seq","a");
  @printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=FWD\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n",2*(index-1)+1,MOVE_Axis,theta_current,2*index,DataFiles_dirgene,RSSFile);
  close(ptf);
  index = index + 1;
end

 delta_phi = delta_phi_RJ;

	j=indmin(abs.(tab_increment_phi-delta_phi));
	ScanInc = tab_increment_phi[j];

	if (facility.=="CAMILL_CATR")||(facility.=="CAMILL_FF")
		ScanStart = -180;
		ScanStop = 180 + ScanStart;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	else
		ScanStart = 0;
		ScanStop = 360 + ScanStart - ScanInc;
		nbdata = nbdata + (ScanStop-ScanStart)/ScanInc + 1;
	end

	DataFiles=string(DataFiles_dirgene,DataFiles_prefixe,@sprintf("_scan%d.amp",index));
	RSSFile = @sprintf("config%d.rss",index);
	ptf = open(RSSFile,"w");
	@printf(ptf,"%s",entete);
	@printf(ptf,"ScanStart= %f\nScanStop= %f\nScanInc= %f\n",ScanStart, ScanStop, ScanInc);
	@printf(ptf,"%s",str_polar);
	@printf(ptf,"%s",str_freq);
	@printf(ptf,"DataFiles= %s\n#End",DataFiles);
	close(ptf);

	ptf = open("SEQ.seq","a");
	@printf(ptf,"[SEQ%d]\nACTION=MOVE\nAXIS=%d\nTARGET=%f\nDIR=REV\n[SEQ%d]\nACTION=RSL\nRSS=%s%s\nCALIBRATED=NO\nCALFILE=\n",2*(index-1)+1,MOVE_Axis,azimuth_offset,2*index,DataFiles_dirgene,RSSFile);
	close(ptf);

	index = index + 1;

writedlm("igloo.txt",igloo);

#/ Transmission measurement estimation
tsweep = 1.0;
tdata_transfer = 0.2;

#// ANT32 config file loading time
tconfig_file = 2;

accelaration_phi = acc_scan;
vitesse_phi = v0_scan;
acceleration_polar = acc_polar;
vitesse_polar = v0_polar;
acceleration_theta = acc_move;
vitesse_theta = v0_move;

tvnom_phi = vitesse_phi/accelaration_phi;

tmp = igloo[:,1];

theta= unique(tmp);
dtheta = abs(theta[2]-theta[1]);
t0 = vitesse_theta/acceleration_theta;
d0 = 0.5*acceleration_theta*t0^2;
tdtheta = 2*t0 + (dtheta-2*d0)/vitesse_theta;

dpolar = 90;
t0 = vitesse_polar/acceleration_polar;
d0 = 0.5*acceleration_polar*t0^2;
tdpolar = 2*t0 + (dpolar-2*d0)/vitesse_polar;

estimated_time = -tdtheta;
for itheta=1:length(theta)
    itmp = find(abs.(igloo[:,1]-theta[itheta]).<0.001);
    phi = igloo[itmp,2];

    dphi = abs(phi[2]-phi[1]);
    nbphi = length(phi)-1;

    t0 = vitesse_phi/accelaration_phi;
    d0 = 0.5*accelaration_phi*t0^2;
    tdphi = 2*t0 + (dphi-2*d0)/vitesse_phi;

    estimated_time = estimated_time + (nbphi*(tdphi+tsweep) + tdpolar + tdata_transfer)*2+tdtheta+tconfig_file;

end


s = mod(estimated_time,60);
m = mod((estimated_time-s)/60,60);
h = (estimated_time-s-60*m)/3600;

@printf("\nEstimated time = %dh %2.2dm %2.2ds\n",h,m,s);

end



#######################################################################
## 			read_sumup_igloo
##
## Read sumup igloo information file
##
##	Input : file_name
##	Output :
##		Facility : facility used for measurement
##		theta_max : maximum theta angular values
##		shift_theta : theta offset
##		shift_phi : phi offset
##		freq_min : min frequency value
##		freq_max : max frequency value
##		nb_freq : nb of frequency values
##		itemp : periodicity of thermal drift estimation procedure
##
#######################################################################



function read_sumup_igloo(file_name="sumup_igloo.txt")

ptf=open(file_name,"r");
data=readdlm(ptf,String);
close(ptf);


Facility=data[find(data.=="Facility"),3][1];
theta_max=floor(1e3*convert_str2float(data[find(data.=="theta_max"),3][1]))*1e-3;
shift_theta=floor(1e3*convert_str2float(data[find(data.=="azimuth_offset"),3][1]))*1e-3;
shift_phi=floor(1e3*convert_str2float(data[find(data.=="roll_offset"),3][1]))*1e-3;
freq_min=floor(1e3*convert_str2float(data[find(data.=="freq_min"),3][1]))*1e-3;
freq_max=floor(1e3*convert_str2float(data[find(data.=="freq_max"),3][1]))*1e-3;
nb_freq=round(Int64,floor(1e3*convert_str2float(data[find(data.=="nb_freq"),3][1]))*1e-3);
itemp=round(Int64,convert_str2float(data[find(data.=="tdrift_periodicity"),3][1]));

tab_freq = collect(linspace(freq_min,freq_max,nb_freq));

return Facility,theta_max, shift_theta, shift_phi,freq_min, freq_max, nb_freq, itemp;

end
