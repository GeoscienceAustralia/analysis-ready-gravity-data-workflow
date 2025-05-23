%==========================================================================
%REFERENCES 
%==========================================================================
%
%Bucha, B., Janak, J., 2014. A MATLAB-based graphical user interface
%program for computing functionals of the geopotential up to ultra-high 
%degrees and orders: Efficient computation at irregular surfaces. Computers 
%and Geosciences 66, 219-227, doi: 10.1016/j.cageo.2014.02.005.
%
%blazej.bucha@stuba.sk, blazej.bucha@gmail.com, juraj.janak@stuba.sk
%
%
%
%==========================================================================
%DOWNLOAD (Version 1.1.2)
%==========================================================================
%
%The source code of the isGrafLab and the pdf file 
%"Definition_of_functionals_of_the_geopotential_used_in_GrafLab_software.pdf"
%are available from:
%http://www.svf.stuba.sk/en/departments/department-of-theoretical-geodesy/science-and-research/downloads.html?page_id=4996
%
%
%
%==========================================================================
%LIST OF CHANGES IN isGrafLab (Irregular Surface GRAvity Field LABoratory)
%==========================================================================
%
%--------------------------------------------------------------------------
%
%May 2014 (isGrafLab 1.00, a modified version of GrafLab 1.1.2)
%
%- Published in Computers & Geosciences 66 (2014), pp. 219-227, 
%  doi: 10.1016/j.cageo.2014.02.005.
%
%--------------------------------------------------------------------------
%
%April 2014 (isGrafLab 1.1)  
%
%- The earlier version of isGrafLab automatically assumed that the spherical
%  harmonic coefficients of the degrees 0 and 1 are set as follows: "C00=1", 
%  "C10=0", "C11=0" and "S11=0". There is no such a restriction in the new 
%  version 1.1. This means that, for example, spherical harmonic coefficients of the 
%  potential of topographic masses can now be imported and used with "nmin = 0" 
%  (instead of "nmin = 2" as in the previous versions of isGrafLab). In general, 
%  in the case of the potential of topographic masses, these coefficietns have 
%  different values from the ones stated above. 
%  Thereby, the following functionals must be computed with "nmin = 0": "Geoid 
%  undulation", "Height anomaly", "Gravity disturbance", and the non-diagonal 
%  elements of the disturbing and the gravitational tensor in the LNOF. The rest 
%  of the functionals can be computed with a value of nmin in the interval 
%  "0 <= nmin <= nmax". For a detailed overview, please, see Table 3 in the source 
%  code "isGrafLab.m". 
%  If the GGM file does not specify any of these coefficients, GrafLab will 
%  automatically use the aforementioned values for the unspecified coefficients. 
%  
%- isGrafLab 1.1 computes the zonal spherical harmonic coefficients of the reference 
%  ellipsoid up to degree "n = 20" (in the previous version this value was set to 
%  n = 10 which is sufficient in the vast majority of practical applications). 
%  The differences compared to "n = 10" are of very small orders of magnitude, 
%  although might be noticeable in the results provided by isGrafLab. For example, 
%  the differences in terms of the disturbing potential are approximately of the 
%  order of 10^-8 m^2*s^-2 or so.
%
%--------------------------------------------------------------------------
%
%August 2015 (isGrafLab 1.1.1)
% 
%- Accelerated routine to compute the fully normalized associated Legendre functions
%  via the extended-range arithmetic approach on 64-bit Matlab. For example, 
%  for "nmax=2190" and 1801 latitudes, the new routine is approximately two times 
%  faster then the previous one. The routine for 32-bit Matlab remains the same.
%
%--------------------------------------------------------------------------  
%
%April 2016 (isGrafLab 1.1.2)
%
%- Minor modifications to ensure compatibility with the latest releases of
%  Matlab.
%
%==========================================================================
%USER MANUAL
%==========================================================================
%
%The GUI of the isGrafLab is visually divided into three panels:
%
%--------------------------------------------------------------------------
%(i) GEOPOTENTIAL MODEL AND REFERENCE SYSTEM SELECTION: At
%first, the input GGM file must be imported using the "Browse. . ." 
%button. The input GGM file must have one of the two standardized 
%structures, see Table 1 and Table 2. In addition to the spherical harmonic
%coefficients, the input file may or may not contain the fifth and
%the sixth column with their standard deviations. isGrafLab is also
%capable to read the input GGM file in a standard format defined by
%ICGEM (International Centre for Global Earth Models). In this case,
%the input file must have the suffix ".gfc" and the spherical harmonic
%coefficients sorted with respect to Table 1 or Table 2. If the input 
%GGM file has this particular structure, GrafLab uses the values of "GM" 
%and "R" from this file, and ignores the values of these two variables
%entered in the GUI. The input file may either be an ASCII file or a 
%binary MAT-file. In case of the GGM with a high maximum degree of SHE, 
%it is recommended to use the binary MAT-file, since it can be loaded 
%much faster.
%
%Table 1: Structure of the input GGM file - spherical harmonic 
%coefficients sorted primarily according to degrees.
%----------------------------------------
%  n   m       C_nm           S_nm
%----------------------------------------
%  2   0   -0.48417E-03    0.00000E+00
%  2   1   -0.20662E-09    0.13844E-08
%  2   2    0.24394E-05   -0.14003E-05
%  3   0    0.95716E-06    0.00000E+00
%----------------------------------------
%
%Table 2: Structure of the input GGM file - spherical harmonic 
%coefficients sorted primarily according to orders.
%----------------------------------------
%  n   m       C_nm           S_nm
%----------------------------------------
%  2   0   -0.48417E-03    0.00000E+00
%  3   0    0.95712E-06    0.00000E+00
%  4   0    0.53998E-06    0.00000E+00
%  5   0    0.68658E-07    0.00000E+00
%----------------------------------------
%
%Most of GGMs have the same values of the geocentric gravitational
%constant and the radius of the reference sphere, therefore in this
%panel isGrafLab automatically offers them for the computation.
%However, they may be simply replaced by the required values, 
%if necessary. Using the arrays "nmin" and "nmax", integer values
%in the intervals nmin 'in' <0,nmax> and nmax 'in' <2,M> may be
%entered (note that there are a few exceptions where the "nmin" value
%is fixed to 0 and cannot be changed, see Table 3 in (iii) "Calculated
%parameters and output selection" below). From the pop-up menu 
%"Ellipsoid", the normal gravity field generated by the 
%equipotential ellipsoid WGS84 (NIMA, 2000) or GRS80 (Moritz, 2000)
%can be selected.
%
%--------------------------------------------------------------------------
%(ii) POINT TYPE SELECTION:  At first, the "Type of the input coordinates"
%(ellipsoidal/spherical) must be specified. The grid must be defined 
%by using the seven self-explanatory arrays in the bottom left of this 
%panel ("Lat." denotes the latitude and "Lon." is the longitude). 
%The array "Height above the reference surface (m)" denotes the constant 
%height of the grid above the reference ellipsoid (GRS80/WGS84) in the 
%case of the ellipsoidal type of the coordinates or above the reference 
%sphere with the radius "R", defined by the GGM, in the case of the 
%spherical coordinates. To this surface refer the points 
%"P_0(r_0^E,theta_0^E,lambda) 'in' E" or "P_0(r_0,theta,lambda) 'in' S" 
%that will be continued to the irregular surface. By a rule of thumb, 
%to this array one can enter the average value computed from the minimum 
%and the maximum height of the irregular surface above the reference surface 
%(the ellipsoid GRS80/WGS84 or the sphere with the radius "R"). The entries 
%must be either in the form of floating point numbers with decimal points 
%or integer values. Latitudes must be entered within the <90�,90�> interval
%and longitudes within the <0�,360�> or <-180�,180�> interval.
%
%Integer value of the Taylor series order ("K >= 0") must be entered 
%into the array "Order of Taylor series:". If "K=0", the functionals 
%are not continued from the regular surface.
%
%Heights of the grid points at the irregular surface must be loaded 
%using the button "Browse...". In case of the ellipsoidal typeof the input 
%coordinates this file must contain ellipsoidal heights. If the spherical 
%coordinates have been chosen, spherical radii of the grid points must be 
%given. In both cases, the input file may have two forms:
%
%a) a matrix
%      -                                              -
%     |(theta_n,lambda_1)   . . .    (theta_n,lambda_m)|
%     |       .             .                .         |
%     |       .               .              .         |
%     |       .                 .            .         |
%     |(theta_1,lambda_1)   . . .    (theta_1,lambda_m)|
%      -                                              -
%
%b) a column vector
%
%      -                -
%     |(theta_1,lambda_1)|
%     |        .         |
%     |        .         |
%     |        .         |
%     |(theta_n,lambda_1)|
%     |        .         |
%     |        .         |
%     |        .         |
%     |(theta_1,lambda_m)|
%     |        .         |
%     |        .         |
%     |        .         |
%     |(theta_n,lambda_m)|
%      -                -
%
%where the point with the coordintates "(theta_1,lambda_1)" is the 
%southernmost and the westernmost point of the grid (arrays 
%"Lat. min (�)" and "Lon. min (�)"), and the point "(theta_n,lambda_m)" 
%is the northernmost and the easternmost point (arrays "Lat. max 
%(�)" and "Lon. max (�)"). The input file may either be an ASCII file 
%or a binary MAT-file. Points with undefined heights/spherical radii 
%should be indicated by value -9999 or NaN. The output values are 
%correspondingly. 
%
%--------------------------------------------------------------------------
%(iii) CALCULATED PARAMETERS AND OUTPUT SELECTION: Using the
%four pop-up menus on the left side of this panel, user can simply
%choose, which functionals of the geopotential are to be computed.
%Note that at least one and maximum four functionals may be computed 
%simultaneously. The summary of the functionals that can be computed in 
%isGrafLab is shown in Table 3. In order to stay brief, we do not introduce 
%here the mathematical formulae for computing each functional, but these 
%can be found in the pdf file
%"Definition_of_functionals_of_the_geopotential_used_in_GrafLab_software.pdf". 
%For evaluating disturbing and gravitational tensor in the LNOF, we
%used the non-singular expressions, which can be found e.g. in Petrovskaya
%and Vershkov (2006). For practical reasons, we slightly modified these
%formulae. The modified formulae can be found in the same pdf file.
%
%Table 3: Functionals of the geopotential available in isGrafLab.
%Explanation of the symbols in the table: "V" - gravitational potential,
%"W" - gravity potential, "g" - gravity, "T" - disturbing potential, 
%"delta g" - gravity disturbance, "DELTA g" - gravity anomaly, "xi" -
%north-south component of deflection of the vertical, "eta" - east-west
%component of deflection of the vertical, "THETA" - total deflection of the
%vertical, "N" - geoid undulation, "zeta_Ell" - generalized height anomaly,
%"zeta" - height anomaly; the subscript "sa" denotes the spherical
%approximation of the functional; ("r", "theta", "lambda") stands for the
%spherical coordinates; ("x","y","z") denotes the coordinates in the local
%north-oriented reference frame; the subscripts "r", "theta", "lambda", 
%"x", "y", "z" and their combinations stand for the derivatives of the 
%functionals with respect to the particular coordinate; the number in the 
%superscript denotes computational demand (computation time of the 
%functional and memory usage during the computation) - (1) small, 
%(2) medium, (3) high, (4) very high; (*) denotes the functionals 
%for which the value of "nmin" cannot be larger than 0.
%-------------------------------------------------------------------------------------------------------------------
%          Actual field                    Disturbing field          Geometrical characteristics of the actual field
%-------------------------------------------------------------------------------------------------------------------
%              V(1)                             T(1)                           xi(2)
%   V_rr(3),V_phiphi(3),V_ll(3)     T_rr(3),T_phiphi(3),T_ll(3)                eta(1)
%   V_rphi(3),V_rl(3),V_phil(3)     T_rphi(3),T_rl(3),T_phil(3)                THETA(2)
%   V_xx(3),V_yy(3),V_zz(3)         T_xx(3),T_yy(3),T_zz(3)                    (*)N(2)
%(*)V_xy(4),(*)V_xz(4),(*)V_yz(4)   (*)T_xy(4),(*)T_xz(4),(*)T_yz(4)           zeta_Ell(1)
%              W(1)                          (*)delta g(3)                     (*)zeta(2)
%              g(2)                          delta g_sa(1)
%              g_sa(1)                       DELTA g_sa(1)
%              W_rr(1)                          T_rr(1)
%-------------------------------------------------------------------------------------------------------------------
%
%To compute geoid undulation "N" and height anomaly "zeta", the
%digital terrain model, e.g. DTM2006.0 (Pavlis et al., 2007),
%must be imported. Only one particular structure of the DTM
%file, shown in Table 1, can be recognized by the isGrafLab.
%If these two functionals are to be computed, immediately after
%clicking the "OK" button, the dialog window from which the input
%DTM file must be imported will appear.
%
%Here, a short note is necessary. The approaches for computing the 
%functionals "Geoid undulation" and "Height anomaly" are the same in 
%isGrafLab as in GrafLab. The reason is that these functionals are referred 
%to particular surfaces, and can not be, by definition, shifted arbitrarily 
%in radial or vertical direction. Moreover, these functionals are computed 
%iteratively as the particular surfaces are not known prior to computation. 
%However, for the convenience of the users we left these functionals also 
%in isGrafLab.
%
%Each functional of the geopotential may be evaluated using
%any of the three approaches for computing fnALFs except for the 
%gravitational and disturbing tensors in the LNOF. Since these 
%non-singular expressions have been slightly modified, the modified 
%forward column method combined with Horner뭩 scheme is not effcient 
%for the new formulae and therefore it was not used in this case.
%
%By clicking the button "Computation of fnALFs", a new dialog
%window will appear, in which user may choose one of the three 
%approaches for evaluating values of fnALFs.
%
%If the Mapping toolbox of MATLAB is installed, computed data may be 
%depicted on a map using automatically selected cartographic projection 
%(e.g. pseudocylindical Robinson projection, equidistant conic projection, 
%equidistant azimuthal projection, ...). By clicking the button "Display 
%data settings", another dialog window will appear. Here, user can set 
%up the required output parameters of the exported map. This option is 
%available only if the computation on a regular grid has been chosen.
%
%The button "Output folder and file" permits to specify the output
%folder and prefix of the all exported files, i.e. without any
%suffix (e.g. "Prefix"). The data file (e.g. "Prefix.txt") with the computed
%data may be created by selecting the checkbox "Export data". 
%The report file, which contains the informations about the
%computation, may be created by selecting the
%"Export report" checkbox. This file automatically obtains name
%with the suffix "_Report.txt", e.g. "Prefix_Report.txt". If the "Display
%data" checkbox has been selected, isGrafLab creates also a
%graphical file (or files, depending on the number of computing
%functionals) according to chosen graphic file format (bmp, emf,
%eps, jpeg, pdf, png or tiff).
%
%When all the required input parameters and input files have
%been entered, after clicking the "OK" button, the computation
%will start. On the left from this button, there is a status line,
%which provides short explanations during the whole computational
%process ("Loading GGM file..." , current value of the variable
%m in the order-dependent loop, "Displaying data..." , etc. ), so that
%a user can clearly see in which part of the computation the isGrafLab is.
%After successful computation, the status "Computation has been
%finished" will appear. If any of the input parameters or input
%files have been entered in a wrong format, isGrafLab will open a
%message dialog or error dialog with description of the error.

function isGrafLabJM(vstpar,GGMname,Latmin,Latmax,Longmin,Longmax,Step,HeightAbove,loadname,F1,F2,F3,F4,outname,outadresar)
R1=0.8; G1=0.8; B1=0.8;
R2=0.95; G2=0.95; B2=0.95;

if vstpar==0  
    
    %Main window 
    M=figure('units','pixels','numbertitle','off','name','isGrafLab 1.1.2',...
        'color',[R1 G1 B1],'position',[300 100 600 600],...
        'tag','okno','menubar','none');      
    a=0.01; b=0.04; c=0.045; d=-0.02;  
    
    %Panels
    %======================================================================
    %Geopotential model and reference system selection panel
	uipanel('Units','normalized','position',[0.06 0.77 0.88 0.21],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'GMaRSSpanel');
    
    %Irregular surface selection panel
	uipanel('Units','normalized','position',[0.06 0.355 0.88 0.4],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'ISSpanel');
    
    %Calculated parameters and output selection panel
	uipanel('Units','normalized','position',[0.06 0.08 0.88 0.26],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'CPaOSpanel');
    
    %Geopotential model and reference system selection
    %======================================================================
    uicontrol('Units','normalized','position',[0.08 0.865+a 0.14 0.035],...
        'style','pushbutton','string','Browse...','tag','GGM',...
        'callback','isGrafLab import_GGM'); %Browse... button
    uicontrol('Units','normalized','position',[0.6 0.865+a 0.32 0.035],...
        'style','checkbox','string','Use maximum degree of GGM',...
        'value',1,'backgroundcolor',[R1 G1 B1],'tag','use'); %Use maximum degree of GGM
    uicontrol('Units','normalized','position',[0.08 0.91+a 0.375 0.025],...
        'style','text','string','Global geopotential model of the Earth',...
        'backgroundcolor',[R1 G1 B1]); %Text Global geopotential model of the Earth
    uicontrol('Units','normalized','position',[0.075 0.82+a 0.22 0.025],...
        'style','text','string','GM of GGM (m3.s-2)',...
        'backgroundcolor',[R1 G1 B1]); %Text GM of GGM (m^3.s^-2)
    uicontrol('Units','normalized','position',[0.08 0.78+a 0.21 0.035],...
        'style','edit','string','3986004.415E+8','backgroundcolor',...
        [R2 G2 B2],'tag','GM'); %Value of GM
    uicontrol('Units','normalized','position',[0.325 0.82+a 0.2 0.025],...
        'style','text','string','R of GGM (m)',...
        'backgroundcolor',[R1 G1 B1],'tag','R_text'); %Text R of GGM (m)
    uicontrol('Units','normalized','position',[0.32 0.78+a 0.21 0.035],...
        'style','edit','string','6378136.3','backgroundcolor',[R2 G2 B2],...
        'tag','R'); %Value of R
    uicontrol('Units','normalized','position',[0.565 0.82+a 0.06 0.025],...
        'style','text','string','nmin','backgroundcolor',[R1 G1 B1],...
        'tag','text_nmin'); %Text nmin
    uicontrol('Units','normalized','position',[0.56 0.78+a 0.08 0.035],...
        'style','edit','string','0','backgroundcolor',[R2 G2 B2],...
        'tag','nmin');  %Value of nmin
    uicontrol('Units','normalized','position',[0.675 0.82+a 0.06 0.025],...
        'style','text','string','nmax','backgroundcolor',[R1 G1 B1],...
        'tag','text_nmax'); %Text nmax
    uicontrol('Units','normalized','position',[0.67 0.78+a 0.08 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','nmax'); %Value of nmax
    uicontrol('Units','normalized','position',[0.795 0.82+a 0.09 0.025],...
        'style','text','string','Ellipsoid','backgroundcolor',[R1 G1 B1],...
        'tag','ell_text'); %Text Ellipsoid
    uicontrol('Units','normalized','position',[0.78 0.717+a 0.13 0.1],...
        'style','popup','string','WGS84|GRS80','backgroundcolor',...
        [R2 G2 B2],'tag','ell','Value',2); %Ellipsoid - pop-up menu
    uicontrol('units','normalized','position',[0.25 0.865+a 0.31 0.035],...
        'style','edit','tag','nameGGM','enable','off','String',GGMname); %Name of the imported GGM file
    uicontrol('Units','normalized','position',[0.06 0.03 0.28 0.025],...
        'style','text','backgroundcolor',[R1 G1 B1],'tag','hlasky'); %Status line

    %Irregular surface selection
    %======================================================================
    
    %Text Input coordinates
    uicontrol('Units','normalized','position',[0.08 0.65+b 0.28 0.025],...
        'style','text','string','Type of the input coordinates:','backgroundcolor',...
        [R1 G1 B1]);    
    
    %Radio button group "Type of the input coordinates"
    c0 = uibuttongroup('visible','on','units','normalized',...
                'Position',[0.375 0.618+b 0.4 0.06],'bordertype','none',...
                'backgroundcolor',[R1 G1 B1],'tag','coordinates'); 
    
    %Radio button: Ellipsoidal
    c1 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.0 0.5 0.4 0.5],'parent',c0,'HandleVisibility','on',...
                'backgroundcolor',[R1 G1 B1],'tag','rbutton1coord'); 
    set(c1,'String','Ellipsoidal','fontname','cambria','fontsize',10);
    
    %Radio button: Spherical
    c2 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.388 0.5 0.4 0.5],'parent',c0,'HandleVisibility','on',...
                'backgroundcolor',[R1 G1 B1],'tag','rbutton2coord'); 
    set(c2,'String','Spherical','fontname','cambria','fontsize',10); 
        
    %Text Order of Taylor series:
    uicontrol('Units','normalized','position',[0.08 0.595+b 0.208 0.025],...
        'style','text','string','Order of Taylor series:','backgroundcolor',...
        [R1 G1 B1]);
    
    %Value Order of Taylor series
    uicontrol('Units','normalized','position',[0.30 0.59+b 0.075 0.035],...
        'style','edit','string','3','backgroundcolor',[R2 G2 B2],...
        'tag','TR');
    
    %Grid
    %----------------------------------------------------------------------
    
    %phi
    disp(num2str(Latmin))
    uicontrol('Units','normalized','position',[0.08 0.5+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fimin','String',num2str(Latmin)); %phi min
    uicontrol('Units','normalized','position',[0.225 0.5+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fistep','String',num2str(Step)); %phi step
    uicontrol('Units','normalized','position',[0.37 0.5+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fimax','String',num2str(Latmax)); %phi max

    %lambda
    uicontrol('Units','normalized','position',[0.08 0.418+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdamin','String',num2str(Longmin)); %lambda min
    uicontrol('Units','normalized','position',[0.225 0.418+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdastep','String',num2str(Step)); %lambda step
    uicontrol('Units','normalized','position',[0.37 0.418+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdamax','String',num2str(Longmax)); %lambda max

    %h
    uicontrol('Units','normalized','position',[0.08 0.335+b 0.4 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','hgrid','String',num2str(HeightAbove)); %h
    
    %Text for latitudes
    uicontrol('Units','normalized','position',[0.075 0.545+b 0.12 0.025],...
        'style','text','string','Lat. min (�)','backgroundcolor',...
        [R1 G1 B1],'tag','fimin_string'); %Text Lat. min (�)
    uicontrol('Units','normalized','position',[0.22 0.545+b 0.12 0.025],...
        'style','text','string','Lat. step (�)','backgroundcolor',...
        [R1 G1 B1],'tag','fistep_string'); %Text Lat. step (�)
    uicontrol('Units','normalized','position',[0.365 0.545+b 0.12 0.025],...
        'style','text','string','Lat. max (�)','backgroundcolor',...
        [R1 G1 B1],'tag','fimax_string'); %Text Lat. max (�)
    
    %Text for longitudes
    uicontrol('Units','normalized','position',[0.075 0.463+b 0.12 0.025],...
        'style','text','string','Lon. min (�)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. min (�)
    uicontrol('Units','normalized','position',[0.22 0.463+b 0.12 0.025],...
        'style','text','string','Lon. step (�)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. step (�)
    uicontrol('Units','normalized','position',[0.365 0.463+b 0.12 0.025],...
        'style','text','string','Lon. max (�)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. max (�)
        
    %Text Height above the reference surface (m)
    uicontrol('Units','normalized','position',[0.08 0.38+b 0.4 0.025],...
        'style','text','string','Height above the reference surface (m)',...
        'backgroundcolor',[R1 G1 B1],'tag','h_string'); 

    %Text: Irregular surface file
    uicontrol('Units','normalized','position',[0.55+d 0.545+b 0.195 0.025],...
        'style','text','string','Irregular surface file','backgroundcolor',[R1 G1 B1],...
        'tag','ell_coord');
    
    %Browse... button
    uicontrol('Units','normalized','position',[0.55+d 0.5+b 0.14 0.035],...
        'style','pushbutton','string','Browse...','tag','import',...
        'callback','isGrafLab input'); 
           
    uicontrol('units','normalized','position',[0.71+d 0.5+b 0.23 0.035],...
        'style','edit','tag','nameDTM','enable','off','String',loadname); %Name of the imported DTM file
    
    %Radio button group "Loading data for spherical or ellipsoidal coordinates"
    %----------------------------------------------------------------------
    
    d0 = uibuttongroup('visible','on','units','normalized',...
                'Position',[0.57+d 0.335+b 0.38 0.08],'bordertype','none',...
                'backgroundcolor',[R1 G1 B1],'tag','Loading_data');
      
    %Text: Structure of the input file:
    uicontrol('Units','normalized','position',[0.55+d 0.425+b 0.24 0.025],...
        'style','text','string','Structure of the input file:','backgroundcolor',[R1 G1 B1],...
        'tag','ell_coord');
               
    %Radio button: Matrix
    d1 = uicontrol('units','normalized','Style','Radio','pos',...
                [0 0.64 0.25 0.3],'parent',d0,'HandleVisibility','on',...
                'backgroundcolor',[R1 G1 B1],'tag','rbutton1Load_data'); 
    set(d1,'String','Matrix','fontname','cambria','fontsize',10,...
        'value',true);
    
    %Radio button: Column vector
    d2 = uicontrol('units','normalized','Style','Radio','pos',...
                [0 0.07 0.45 0.3],'parent',d0,'HandleVisibility','on',...
                'backgroundcolor',[R1 G1 B1],'tag','rbutton2Load_data'); 
    set(d2,'String','Column vector','fontname','cambria','fontsize',10,...
        'value',false);
       
    %Calculated parameters and output selection
    %======================================================================
    
    %The first functional
    uicontrol('Units','normalized','position',[0.08 0.16+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P1','Value',F1);
    %The second functional
	uicontrol('Units','normalized','position',[0.08 0.105+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P2','Value',F2);    
    %The third functional
	uicontrol('Units','normalized','position',[0.08 0.05+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P3','Value',F3);   
    %The fourth functional
	uicontrol('Units','normalized','position',[0.08 -0.005+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P4','Value',F4);

    %----------------------------------------------------------------------
     
    %Computation of fnALFs
    uicontrol('Units','normalized','position',[0.42 0.222+c 0.25 0.038],...
        'style','pushbutton','string','Computation of fnALFs',...
        'callback','isGrafLab fnALFs','tag','fnALFs'); 
    
    %Display data settings
    uicontrol('Units','normalized','position',[0.42 0.167+c 0.25 0.038],...
        'style','pushbutton','string','Display data settings',...
        'callback','isGrafLab Display_data_settings','tag','DDS');
    
    %Output folder and file
    uicontrol('Units','normalized','position',[0.42 0.112+c 0.25 0.038],...
        'style','pushbutton','string','Output folder and file',...
        'callback','isGrafLab Output_folder','tag','outfolder');
        
    %Export data checkbox
    uicontrol('Units','normalized','position',[0.72+d 0.228+c 0.2 0.025],...
        'style','checkbox','string','Export data','tag','export',...
        'callback','isGrafLab output','backgroundcolor',[R1 G1 B1],'value',0);
    
    %Export report checkbox
    uicontrol('Units','normalized','position',[0.72+d 0.175+c 0.2 0.025],...
        'style','checkbox','string','Export report','backgroundcolor',...
        [R1 G1 B1],'tag','report','value',0);
    
    %Export data in *.mat
    uicontrol('Units','normalized','position',[0.72+d 0.12+c 0.23 0.025],...
        'style','checkbox','string','Export data in *.mat','backgroundcolor',...
        [R1 G1 B1],'tag','datamat','value',1);
    
    %OK button
    uicontrol('Units','normalized','position',[0.35 0.015 0.13 0.05],...
        'style','pushbutton','string','OK','Callback','isGrafLab OK');
    
    %Close button
	uicontrol('Units','normalized','position',[0.55 0.015 0.13 0.05],...
        'style','pushbutton','string','Close','Callback','isGrafLab Close');

    %Set font to Cambria
    set(get(M,'children'),'fontname','cambria','fontsize',10);  
    set(findobj('tag','ISSpanel'),'title','Irregular surface selection',...
        'fontname','cambria','fontsize',8);   
    set(findobj('tag','GMaRSSpanel'),'title',...
        'Geopotential model and reference system selection','fontname',...
        'cambria','fontsize',8); 
    set(findobj('tag','CPaOSpanel'),'title',...
        'Calculated parameters and output selection','fontname',...
        'cambria','fontsize',8);
    set(findobj('tag','R_text'),'userdata',outname);
    set(findobj('tag','ell_text'),'userdata',outadresar);
    set(findobj('tag','R'),'userdata',GGMname);
    set(findobj('tag','ell'),'userdata','');
    set(findobj('tag','use'),'userdata',loadname);
    set(findobj('tag','datamat'),'userdata','');    
    isGrafLabJM OK
else    
    switch(vstpar)
        
        case('import_GGM') %Click on the Browse... button in the Geopotential model and reference system selection panel
            
                [GGMname,GGMadresar]=uigetfile('*.*','Select GGM File');
                if GGMname==0
                else
                    set(findobj('tag','R'),'userdata',GGMname);
                    set(findobj('tag','ell'),'userdata',GGMadresar);
                    
                    set(findobj('tag','nameGGM'),'string',GGMname); %Display the name of the imported GGM file
                end
                
        case('input') %Click on the Browse... button in the Irregular surface panel
            
            [loadname,loadadresar]=uigetfile('*.*','Select File Containing Heights of Irregular Surface');
            if loadname==0
            else
                set(findobj('tag','use'),'userdata',loadname);
                set(findobj('tag','datamat'),'userdata',loadadresar);
       
                set(findobj('tag','nameDTM'),'string',loadname); %Display the name of the imported GGM file
            end
        
%         case('Output_folder') %Click on the Output folder and file button
%             
%             [outname,outadresar]=uiputfile('*.*');
%             if outname==0
%             else
%                 if find(outname=='.')>0
%                     outname=outname(1:(find(outname=='.')-1));
%                 end
%                 set(findobj('tag','R_text'),'userdata',outname);
%                 set(findobj('tag','ell_text'),'userdata',outadresar);
%             end
            
        case('Display_data_settings') %Click on the Display data settings
            
            %Main window
            D=figure('units','pixels','numbertitle','off','name',...
                'Display data settings','color',[R1 G1 B1],'position',...
                [320 150 550 350],'tag','oknoDDS','menubar','none');
            
            %Panel
            uipanel('Units','normalized','position',[0 0 0.26 1],...
                'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],...
                'fontsize',8,'fontname','cambria');
            
            %Display data checkbox
            uicontrol('Units','normalized','position',[0.05 0.87 0.2 0.05],...
                'style','checkbox','string','Display data','backgroundcolor',...
                [R1 G1 B1],'tag','Display');
            
            display_data=get(findobj('tag','DDS'),'userdata');
            if display_data==0
            elseif display_data==1
                set(findobj('tag','Display'),'value',1);
            end

            %Text next to the Display data checkbox
            g0=uicontrol('Units','normalized','position',[0.3 0.83 0.65 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'To export a graphic file, select this checkbox. The data will be depicted on a map using automatically selected cartographic projection.');
            set(g0,'HorizontalAlignment','left');
            
            %Text Graphic format
            uicontrol('Units','normalized','position',[0.05 0.75 0.16 0.05],...
                'style','text','string','Graphic format','backgroundcolor',...
                [R1 G1 B1],'tag','format');
            
            %Graphic format pop-up menu
            uicontrol('Units','normalized','position',...
                [0.05 0.64 0.16 0.1],'style','popup','string',...
                '*.bmp|*.emf|*.eps|*.jpeg|*.pdf|*.png|*.tiff',...
                'backgroundcolor',[R2 G2 B2],'tag','pripona');
            
            prip=get(findobj('tag','nmin'),'userdata');     
            if isempty(prip)
                set(findobj('tag','pripona'),'value',6);
            else  
                set(findobj('tag','pripona'),'value',prip);
            end
            
            %Text next to the Graphic format file
            g1=uicontrol('Units','normalized','position',[0.3 0.67 0.65 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Select one of the graphic format files. For a vector output it is recommended to use *.eps graphic file and *.png format for a bitmap output.');
            set(g1,'HorizontalAlignment','left');
            
            %Text Colormap
            uicontrol('Units','normalized','position',[0.05 0.57 0.16 0.05],...
                'style','text','string','Colormap',...
                'backgroundcolor',[R1 G1 B1]);
            
            %Colormap pop-up menu
            uicontrol('Units','normalized','position',[0.05 0.46 0.16 0.1],...
                'style','popup','string',...
                'jet|HSV|hot|cool|spring|summer|autumn|winter|gray|bone|copper|pink|lines',...
                'backgroundcolor',[R2 G2 B2],'tag','colormap');
            set(findobj('tag','colormap'),'value',2);
            
            color=get(findobj('tag','nmax'),'userdata');     
            if isempty(color)
                set(findobj('tag','colormap'),'value',1);
            else  
                set(findobj('tag','colormap'),'value',color);
            end
            
            %Text next to the colormap
            g1=uicontrol('Units','normalized','position',[0.3 0.417 0.65 0.2],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Select a colormap of the output file. Mostly it is recommended to use the jet colormap, which ranges from blue to red, and passes through the colors cyan, yellow, and orange.');
            set(g1,'HorizontalAlignment','left');
            
            %Text Number of colors
            uicontrol('Units','normalized','position',[0.04 0.39 0.19 0.05],...
                'style','text','string','Number of colors',...
                'backgroundcolor',[R1 G1 B1]);
            
            %Value of number of colors
            uicontrol('Units','normalized','position',[0.08 0.32 0.1 0.06],...
                'style','edit','string',15,'backgroundcolor',...
                [R2 G2 B2],'tag','skala');
            
            %Text next to the number of colors
            g1=uicontrol('Units','normalized','position',[0.3 0.31 0.65 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1]);
            set(g1,'HorizontalAlignment','left','string',...
                'Enter a number of colors of the selected colormap. Note that processing time may increase to a several minutes, if a large number of colors has been entered for a large data set.');
            
            ncolor=get(findobj('tag','text_nmin'),'userdata');
            if isempty(ncolor)                
            else
                set(findobj('tag','skala'),'string',ncolor);
            end
            
            %Text DPI 
            uicontrol('Units','normalized','position',[0.08 0.24 0.1 0.05],...
                'style','text','string','DPI','backgroundcolor',[R1 G1 B1]);
            
            %Value of DPI
            uicontrol('Units','normalized','position',[0.08 0.17 0.1 0.06],...
                'style','edit','string',300,'backgroundcolor',[R2 G2 B2],...
                'tag','DPI');
            
            DPI=get(findobj('tag','text_nmax'),'userdata');
            if isempty(DPI)                
            else
                set(findobj('tag','DPI'),'string',DPI);
            end
            
            %Text next to the DPI
            g1=uicontrol('Units','normalized','position',[0.3 0.03 0.65 0.2],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Enter a value of dots per inch of the output file.');
            set(g1,'HorizontalAlignment','left');
            
            %OK button
            uicontrol('Units','normalized','position',[0.35 0.02 0.13 0.08],...
                'style','pushbutton','string','OK','Callback','isGrafLab OKDDS');
    
            %Close button
            uicontrol('Units','normalized','position',[0.55 0.02 0.13 0.08],...
                'style','pushbutton','string','Close','Callback',...
                'isGrafLab CloseDDS');
            
            %setting font to Cambria
            set(get(D,'children'),'fontname','cambria','fontsize',10)
            
        case('fnALFs') %Click on the Computation of fnALFs
            
            %Main window
            F=figure('units','pixels','numbertitle','off','name',...
                'Approaches to compute fully normalized associated Legendre functions','color',...
                [R1 G1 B1],'position',[320 150 550 350],'tag',...
                'oknofnALFs','menubar','none'); 
            
            %Radio button group
            u0 = uibuttongroup('visible','on','units','normalized',...
                'Position',[0 0 .26 1],'tag','volbaALFs'); 

            %The first radiobutton
            u1 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.7 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton1'); 
            set(u1,'String','<html>Standard forward<br>column method',...
                'fontname','cambria','fontsize',10);
            note1=uicontrol('Units','normalized','position',...
                [0.3 0.74 0.65 0.2],'style','text','backgroundcolor',...
                [R1 G1 B1]);
            set(note1,'HorizontalAlignment','left','string',...
                'It is recommended to use the standard forward column method for all latitudes up to the maximum degree 1800. However, this method may also be used for the latitudes <0�,56�> and <80�,90�> up to the maximum degree 2190.');
            
            %The second radiobutton
            u2 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.42 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton2'); 
            set(u2,'String',...
                '<html>Modified forward<br>column method<br>combined with<br>Horner큦 scheme',...
                'fontname','cambria','fontsize',10);
            note2=uicontrol('Units','normalized','position',[0.3 0.38 0.65 0.3],...
                'style','text','backgroundcolor',[R1 G1 B1]);
            set(note2,'HorizontalAlignment','left','string',...
                'It is recommended to use the modified forward column method combined with Horner큦 scheme for all latitudes and maximum degrees ranging from 1801 to 2700. This method may also be used for lower degrees than 1801, but cannot be applied to degrees higher than 2700 due to the overflow problem.');            
            
            %The third radiobutton
            u3 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.12 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton3'); 
            set(u3,'String','<html>Extended-range<br>arithmetic',...
                'fontname','cambria','fontsize',10);
            note3=uicontrol('Units','normalized','position',...
                [0.3 0.018 0.65 0.3],'style','text','backgroundcolor',...
                [R1 G1 B1]);
            set(note3,'HorizontalAlignment','left','string',...
                'The extended-range arithmetic approach may be used for all latitudes up to an arbitrary degree essentially.');                        
            
            %OK button
            uicontrol('Units','normalized','position',[0.35 0.02 0.13 0.08],...
                'style','pushbutton','string','OK','Callback','isGrafLab OKfnALFs');
    
            %Close button
            uicontrol('Units','normalized','position',[0.55 0.02 0.13 0.08],...
                'style','pushbutton','string','Close','Callback','isGrafLab ClosefnALFs');
            
            %The chosen approach for computation of fnALFs
            volbaALFs=get(findobj('tag','fnALFs'),'userdata');
            if isempty(volbaALFs)
                volbaALFs=1;
            end
            
            %Mark the chosen radiobutton
            if volbaALFs==1
                set(u1,'value',1);
            elseif volbaALFs==2
                set(u2,'value',1);
            elseif volbaALFs==3
                set(u3,'value',1);
            end
            
            set(get(F,'children'),'fontname','cambria','fontsize',10);
            
        case('ClosefnALFs') %Click on the Close button in the Computation of fnALFs window
            
            close
            
        case('OKfnALFs') %Click on the OK button in the Computation of fnALFs window
            
            volbaALFs=get(findobj('tag','volbaALFs'),'selectedobject');

            if get(findobj('tag','rbutton1'),'value')==1
                volbaALFs=1;
            elseif get(findobj('tag','rbutton2'),'value')==1
                volbaALFs=2;
            elseif get(findobj('tag','rbutton3'),'value')==1
                volbaALFs=3;
            end
            
            set(findobj('tag','fnALFs'),'userdata',volbaALFs);            
            
            close
            
        case('OKDDS') %Click on the OK button in the Display data settings window
            
            display_data=get(findobj('tag','Display'),'value');
            set(findobj('tag','DDS'),'userdata',display_data);  
            
            prip=get(findobj('tag','pripona'),'value');
            set(findobj('tag','nmin'),'userdata',prip);
            
            color=get(findobj('tag','colormap'),'value');
            set(findobj('tag','nmax'),'userdata',color);

            ncolor=get(findobj('tag','skala'),'string');
            ncolor=str2double(ncolor);
            set(findobj('tag','text_nmin'),'userdata',ncolor);
            
            DPI=get(findobj('tag','DPI'),'string');
            DPI=str2double(DPI);
            set(findobj('tag','text_nmax'),'userdata',DPI);
            
            %Check of the entered value of DPI and number of colors
            if isnan(DPI)==1 || DPI<0
                errordlg('Entered DPI value must be larger than zero.',...
                    'Display data settings');
                error('Entered DPI value must be larger than zero.');
            end
            
            if isnan(ncolor)==1
                errordlg('Entered value of number of colors is not correct.',...
                    'Display data settings');
                error('Entered value of number of colors is not correct.');
            end
            
            if rem(DPI,1)~=0
                errordlg('Value of DPI must be a positive integer.',...
                    'Display data settings');
                error('Value of DPI must be a positive integer.');
            end
            
            if ncolor<2
                errordlg('Value of number of colors must be larger than 1.',...
                    'Display data settings');
                error('Value of number of colors must be larger than 1.');
            end
            
            if rem(ncolor,1)~=0
                errordlg('Value of number of colors must be a positive integer.',...
                    'Display data settings');
                error('Value of number of colors must be a positive integer.');
            end
            
            close
            
        case('CloseDDS') %Click on the Close button in the Display data settings window
            
            close
            
        case('OK') %Click on the OK button in the main isGrafLab window 
        
            set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow; 
                      
            %Order of Taylor series
            TR=str2double(get(findobj('tag','TR'),'string'));
            if rem(TR,1)~=0
                errordlg('Order of Taylor series must be a non-negative integer.',...
                    'Irregular surface selection')
                error('Order of Taylor series must be a non-negative integer.')
            end
            if TR<0
                errordlg('Order of Taylor series must be an integer value higher or equal to zero.',...
                    'Irregular surface selection')
                error('Order of Taylor series must be an integer value higher or equal to zero.')
            end
 
            %Ellipsoid
%             if get(findobj('tag','ell'),'value') == 1; %GRS80              
                GMEl=3986005*10^8; %Geocentric gravitational constant of GRS80
                aEl=6378137; %Semimajor axis of GRS80
                eEl=sqrt(0.006694380022903416); %First eccentricity of GRS80
                omegaEl=7292115*10^-11; %Angular velocity of GRS80
                CEl_20=-108263*10^-8/sqrt(5); %Fully normalized C_20 of GRS80
%             elseif get(findobj('tag','ell'),'value') == 2; %WGS84
%                 GMEl=3986004.418*10^8; %Geocentric gravitational constant of WGS84
%                 aEl=6378137; %Semimajor axis of WGS84
%                 fEl=1/298.257223563; %Flattening of WGS84
%                 omegaEl=7292115*10^-11; %Angular velocity of WGS84
%                 CEl_20=-0.484166774985*10^-3; %Fully normalized C_20 of WGS84
%                 eEl=sqrt(fEl*(2-fEl)); %First eccentricity of WGS84
%             end         
            
            %GM and R of GGM
            GM=str2double(get(findobj('tag','GM'),'string'));            
            R=str2double(get(findobj('tag','R'),'string'));
            if GM<=0
                errordlg('Value of GM must be larger than zero.',...
                    'Error in geopotential model and reference system selection')
                error('Value of GM must be larger than zero.')
            end
            if isnan(GM)
                errordlg('Please input the value of GM.',...
                    'Error in geopotential model and reference system selection')
                error('Please input the value of GM.')
            end
            if R<=0
                errordlg('R must be larger than zero.',...
                    'Error in geopotential model and reference system selection')
                error('R value must be larger than zero.')
            end     
            if isnan(R)
                errordlg('Please input the value of R.',...
                    'Error in geopotential model and reference system selection')
                error('Please input the value of R.')
            end
                                           
            %Selection of functionals of the geopotential
            volbapar1=get(findobj('tag','P1'),'value'); %ID of the first functional
            volbapar2=get(findobj('tag','P2'),'value'); %ID of the second functional
            volbapar3=get(findobj('tag','P3'),'value'); %ID of the third functional
            volbapar4=get(findobj('tag','P4'),'value'); %ID of the fourth functional
            volbapar=[volbapar1;volbapar2;volbapar3;volbapar4];
            pocetpar=length(find(volbapar>1));

            %Error messages for selection of functionals of the geopotential
            if pocetpar==1
                if volbapar1>1 && volbapar2==1 && volbapar3==1 && volbapar4==1
                else
                    errordlg('Please select the functional of the geopotential in the first pop-up menu.',...
                        'Calculated parameters and output selection');
                    error('Please select the functional of the geopotential in the first pop-up menu.');
                end
            elseif pocetpar==2
                if volbapar1>1 && volbapar2>1 && volbapar3==1 && volbapar4==1
                    if volbapar1 == volbapar2
                        errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                            'Calculated parameters and output selection');
                        error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                    end
                else
                    errordlg('Please select the functionals of the geopotential in the first and the second pop-up menu.',...
                        'Calculated parameters and output selection');
                    error('Please select the functionals of the geopotential in the first and the second pop-up menu.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 %If geoid or height anomaly is to be computed
                    errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                        'Calculated parameters and output selection');
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            elseif pocetpar==3
                if volbapar1>1 && volbapar2>1 && volbapar3>1 && volbapar4==1
                    if (volbapar1 == volbapar2) || (volbapar1 == volbapar3) || (volbapar2 == volbapar3)
                        errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                            'Calculated parameters and output selection');
                        error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                    end
                else
                    errordlg('Please select the functionals of the geopotential in the first, second and third pop-up menu.',...
                        'Calculated parameters and output selection');
                    error('Please select the functionals of the geopotential in the first, second and third pop-up menu.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 || length(nonzeros(volbapar==10 | volbapar==23))==2 %If geoid or height anomaly is to be computed
                    errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                        'Calculated parameters and output selection');
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            elseif pocetpar==4
                if (volbapar1 == volbapar2) || (volbapar1 == volbapar3) || (volbapar1 == volbapar4) || (volbapar2 == volbapar3) || (volbapar2 == volbapar4) || (volbapar3 == volbapar4)
                    errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                        'Calculated parameters and output selection');
                    error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 || length(nonzeros(volbapar==10 | volbapar==23))==2 || length(nonzeros(volbapar==10 | volbapar==23))==3 %If geoid or height anomaly is to be computed
                    errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                        'Calculated parameters and output selection');
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            end
            
            %Identification of type of the input coordinates
            coord=get(findobj('tag','rbutton2coord'),'value');

            if coord==1 %Entered spherical coordinates                                                
                if any(volbapar==10) || any(volbapar==23) 
                    errordlg('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. Grid that has been entered in the spherical coordinates (constant value of the radius r) does not refer to this surface. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates and set value in the array Height above the reference surface (m) to zero.','Calculated parameters and output selection')
                    error('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. Grid that has been entered in the spherical coordinates (constant value of the radius r) does not refer to this surface. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates and set value in the array Height above the reference surface (m) to zero.')
                end
            end
            
            if volbapar1==1 && volbapar2==1 && volbapar3==1 && volbapar4==1
                errordlg('Please choose at least one functional of the geopotential.',...
                    'Error in calculated paramters and output selection')
                error('Please choose at least one functional of the geopotential.')   
            end 
            
            %Identification of computation of fnALFs approach
            volbaALFs=get(findobj('tag','fnALFs'),'userdata');
            if isempty(volbaALFs)
                volbaALFs=1;
            end
            
            %Identification of Display data
            display_data=get(findobj('tag','DDS'),'userdata');
            if isempty(display_data)
                display_data=0;
            end
            
            %Identification of the output folder and file
            outname=get(findobj('tag','R_text'),'userdata');
            outadresar=get(findobj('tag','ell_text'),'userdata');
            if isempty(outname)
                set(findobj('tag','hlasky'),'string','Select output folder and output file',...
                       'fontsize',8,'foregroundcolor','k'); drawnow;
                
                warn2=warndlg('Output folder and output file were not specified. Click OK and then select an output folder and output file. After the selection, the computation will start.');
                waitfor(warn2);
                
                [outname,outadresar]=uiputfile('*.*');
                if outname==0
                    errordlg('Output folder and output file must be specified!');
                    
                    set(findobj('tag','hlasky'),'string','',...
                       'fontsize',8,'foregroundcolor','k'); drawnow;
                   
                    error('Output folder and output file must be specified!');
                else
                    if find(outname=='.')>0
                        outname=outname(1:(find(outname=='.')-1));
                    end
                end
                
                set(findobj('tag','hlasky'),'string','',...
                       'fontsize',8,'foregroundcolor','k'); drawnow;
            end
            
            %Identification of loaded data
            load_matrix=get(findobj('tag','rbutton1Load_data'),'value');
            load_vector=get(findobj('tag','rbutton2Load_data'),'value');

            %Check of the entered values of DPI a number of colors
            if display_data==1
                ncolor=get(findobj('tag','text_nmin'),'userdata');
                DPI=get(findobj('tag','text_nmax'),'userdata');
                if isnan(DPI)==1 || DPI<0
                    errordlg('Entered DPI value must be larger than zero.',...
                        'Display data settings');
                    error('Entered DPI value must be larger than zero.');
                elseif isnan(ncolor)==1 || ncolor<2
                    errordlg('Value of number of colors must be larger than 1.',....
                        'Display data settings');
                    error('Value of number of colors must be larger than 1.');
                end
            end
            
            %Check, if the computation of tensors in the LNOF using the MFCM 
            %has been selected (not allowed)
            if any(volbapar==8) || any(volbapar==9) || any(volbapar==14) || any(volbapar==15)
                if volbaALFs==2
                    errordlg('The following functionals of the geopotential cannot be computed on a grid using the modified forward column method combined with Horner큦 scheme: Disturbing_tensor_Txx_Tyy_Tzz, Disturbing_tensor_Txy_Txz_Tyz, Gravitational_tensor_Vxx_Vyy_Vzz and Gravitational_tensor_Vxy_Vxz_Vyz. In the case of a high maximum degree (1800 and higher), it is recommended to use the extended-range arithmetic approach. However, for the point-wise computation of the above mentioned functionals, the modified forward column method combined with Horner큦 scheme can be applied.',...
                        'Calculated parameters and output selection');
                    error('The following functionals of the geopotential cannot be computed on a grid using the modified forward column method combined with Horner큦 scheme: Disturbing_tensor_Txx_Tyy_Tzz, Disturbing_tensor_Txy_Txz_Tyz, Gravitational_tensor_Vxx_Vyy_Vzz and Gravitational_tensor_Vxy_Vxz_Vyz. In the case of a high maximum degree (1800 and higher), it is recommended to use the extended-range arithmetic approach. However, for the point-wise computation of the above mentioned functionals, the modified forward column method combined with Horner큦 scheme can be applied.');
                end
            end
                                              
            tic %Start clock to measure computation time
            
            %Loading DMR
            if get(findobj('tag','P1'),'value')==10 || get(findobj('tag','P2'),'value')==10 || get(findobj('tag','P3'),'value')==10 || get(findobj('tag','P4'),'value')==10 || get(findobj('tag','P1'),'value')==23 || get(findobj('tag','P2'),'value')==23 || get(findobj('tag','P3'),'value')==23 || get(findobj('tag','P4'),'value')==23
               set(findobj('tag','hlasky'),'string','Please select DTM file',...
                   'fontsize',8,'foregroundcolor','k'); drawnow;
               
               [loadnameDMR,loadadresarDMR]=uigetfile('*.*','Select File Containg Spherical Harmonic Coefficients of the Topography');

               if loadnameDMR==0
                   errordlg('To compute Geoid_undulation/Height_anomaly, DTM file must be imported!',...
                       'Error in geopotential model and reference system selection');
                   error('To compute Geoid_undulation/Height_anomaly, DTM file must be imported!')
               else
                   set(findobj('tag','hlasky'),'string',...
                       'Loading DTM file...','fontsize',8,...
                       'foregroundcolor','k'); drawnow;

                   if strcmp(loadnameDMR(end-3:end),'.mat')
                       DMR=load([loadadresarDMR,loadnameDMR]);
                       DMR=struct2cell(DMR);
                       DMR=cell2mat(DMR);
                   else
                       DMR=load([loadadresarDMR,loadnameDMR]);
                   end

                   [rows_DMR,cols_DMR]=size(DMR); %#ok<*ASGLU>
                   if cols_DMR<4
                       errordlg('Wrong format of the input DTM file.',...
                           'Geopotential model and reference system selection');
                       error('Wrong format of the input DTM file.')
                   end
                   
                   DMR=DMR(:,1:4);
                   
                   set(findobj('tag','hlasky'),'string','',...
                       'foregroundcolor','k'); drawnow;
               end
            end                                      

            %Loading GGM                              
            GGMname=get(findobj('tag','R'),'userdata');
            GGMadresar=get(findobj('tag','ell'),'userdata');

            if isempty(GGMname) %Error message, if GGM file has not been imported
                errordlg('Please input geopotential model file.',...
                    'Error in global geopotential model and reference system selection');
                error('Please input geopotential model file.')
            end

            set(findobj('tag','hlasky'),'string',...
                    'Loading GGM file...','fontsize',8,...
                    'foregroundcolor','k'); drawnow;

            %Loading GGM
            if strcmp(GGMname(end-3:end),'.gfc') %Input data in ICGEM format 
                fGGMid=fopen([GGMadresar,GGMname]);

                while(~feof(fGGMid))                       
                    s=fgetl(fGGMid); 
                    if strncmpi(s,'earth_gravity_constant',22)
                        GM=str2num(s(23:end)); %#ok<*ST2NM>
                    end
                    if strncmpi(s,'radius',6)
                        R=str2num(s(7:end));
                    end
                    if strncmpi(s,'end_of_head',11)
                        break
                    end
                end

                GGM=textscan(fGGMid,'%s%f%f%f%f%f%f');

                fclose(fGGMid);

                GGM=GGM(2:end);
                GGM=cell2mat(GGM);
            elseif strcmp(GGMname(end-3:end),'.mat') %Input data in MAT format 
                GGM=load([GGMadresar,GGMname]);
                GGM=struct2cell(GGM);
                GGM=cell2mat(GGM);
            else
                GGM=load([GGMadresar,GGMname]);
            end
            
            [rows_GGM,cols_GGM]=size(GGM);
            if cols_GGM<4
                errordlg('Wrong format of the input GGM file.',...
                    'Geopotential model and reference system selection');
                error('Wrong format of the input GGM file.')
            end
            
            GGM=GGM(:,1:4);
            GGM=sortrows(GGM,1);
            stupen=GGM(:,1);
            rad=GGM(:,2);
            C=GGM(:,3);
            S=GGM(:,4);
            nmaxGGM=max(GGM(:,1));

            del00=find(stupen==0 & rad==0,1);
            if isempty(del00) %If the coefficients of the degree 0 and order 0 are missing
                stupen=[0;stupen];
                rad=[0;rad];
                C=[1;C];
                S=[0;S];
            end
                
            del10=find(stupen==1 & rad==0,1); %If the coefficients of the degree 1 and order 0 are missing
            if isempty(del10)
                stupen=[stupen(1);1;stupen(2:end)];
                rad=[rad(1);0;rad(2:end)];
                C=[C(1);0;C(2:end)];
                S=[S(1);0;S(2:end)];
            end                             

            del11=find(stupen==1 & rad==1,1); %If the coefficients of the degree 1 and order 1 are missing               
            if isempty(del11)
                stupen=[stupen(1:2);1;stupen(3:end)];
                rad=[rad(1:2);1;rad(3:end)];
                C=[C(1:2);0;C(3:end)];
                S=[S(1:2);0;S(3:end)];
            end
            clear del00 del10 del11

            %Identification of GGM file format
            if stupen(1)==0 && stupen(2)==1 && stupen(3)==1 && stupen(4)==2 && stupen(5)==2 && rad(1)==0 && rad(2)==0 && rad(3)==1 && rad(4)==0 && rad(5)==1
            else
                errordlg('Wrong format of the input GGM file.',...
                    'Geopotential model and reference system selection');
                error('Wrong format of the input GGM file.')
            end               

            set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;
            
            %Value of nmin and error messages
            nmin=str2double(get(findobj('tag','nmin'),'string'));
            if nmin<0 
                errordlg('Value of nmin cannot be negative.',...
                    'Error in geopotential model and reference system selection')
                error('Value of nmin cannot be negative.')
            elseif nmin>nmaxGGM
                errordlg('Value of nmin exceedes nmax value of GGM.',...
                    'Error in geopotential model and reference system selection')
                error('Value of nmin exceedes nmax value of GGM.')
            end
            if isnan(nmin)==1
                    errordlg('Please input the nmin value.',...
                        'Error in geopotential model and reference system selection')
                    error('Please input the nmin value.')
            end
            if rem(nmin,1)~=0
                errordlg('Value of nmin must be an integer.',...
                    'Error in geopotential model and reference system selection')
                error('Value of nmin must be an integer.')
            end
            
            %Value of nmax and error messages   
            if get(findobj('tag','use'),'value')==1
                nmax=nmaxGGM;
            else
                nmax=str2double(get(findobj('tag','nmax'),'string'));

                if nmax>nmaxGGM
                    errordlg('Entered value of nmax exceedes nmax value of GGM.',...
                        'Error in geopotential model and reference system selection')
                    error('Entered value of nmax exceedes nmax value of GGM.')
                elseif nmax<2
                    errordlg('Value of nmax must be at least 2.',...
                        'Error in geopotential model and reference system selection')
                    error('Value of nmax must be at least 2.')
                elseif nmin>nmax
                    errordlg('Value of nmin cannot be larger than nmax value.',...
                        'Error in geopotential model and reference system selection')
                    error('Value of nmin cannot be larger than nmax value.')
                end
                if isnan(nmax)==1
                    errordlg('Please input the nmax value.',...
                        'Error in geopotential model and reference system selection')
                    error('Please input the nmax value.')
                end
                
                if rem(nmax,1)~=0
                    errordlg('Value of nmax must be an integer.',...
                        'Error in geopotential model and reference system selection')
                    error('Value of nmax must be an integer.')
                end
            end  
            
            if nmin>0
                if any(volbapar==20)
                    errordlg('Gravity disturbance cannot be computed if nmin>0.',...
                        'Calculated parameters and output selection');
                    error('Gravity disturbance cannot be computed if nmin>0.');
                end
            end
                        
            if nmin>0
                if any(volbapar==10) || any(volbapar==23)
                    errordlg('Geoid_undulation and Height_anomaly cannot be computed if nmin>0.');
                    error('Geoid_undulation and Height_anomaly cannot be computed if nmin>0.');
                end
                    
                if any(volbapar==9) || any(volbapar==15)
                    errordlg('The following functionals of the geopotential cannot be computed if nmin>0: Disturbing_tensor_Txy_Txz_Tyz and Gravitational_tensor_Vxy_Vxz_Vyz.',...
                        'Calculated parameters and output selection');
                    error('The following functionals of the geopotential cannot be computed if nmin>0: Disturbing_tensor_Txy_Txz_Tyz and Gravitational_tensor_Vxy_Vxz_Vyz.');
                end
                    
                if any(volbapar==8) || any(volbapar==14) 
                    LNOFnmin=1; %Logical 1
                     
                    if any(volbapar~=8 & volbapar~=9 & volbapar~=14 & volbapar~=15 & volbapar~=1)                                              
                        errordlg('The following functionals of the geopotential cannot be computed simultaneously with other functionals if nmin>0: Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz.',...
                        'Calculated parameters and output selection');
                        error('The following functionals of the geopotential cannot be computed simultaneously with other functionals if nmin>0: Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz.');
                    end
                else
                    LNOFnmin=0; %Logical 0
                end
            end                              
                                                          
            %Sorting spherical harmonic coeffcients of DTM
            if get(findobj('tag','P1'),'value')==10 || get(findobj('tag','P2'),'value')==10 || get(findobj('tag','P3'),'value')==10 || get(findobj('tag','P4'),'value')==10 || get(findobj('tag','P1'),'value')==23 || get(findobj('tag','P2'),'value')==23 || get(findobj('tag','P3'),'value')==23 || get(findobj('tag','P4'),'value')==23
                DMR(DMR(:,1)>nmaxGGM,:)=[];

                if DMR(1,1)==0 && DMR(2,1)==1 && DMR(3,1)==1 && DMR(4,1)==2 && DMR(1,2)==0 && DMR(2,2)==0 && DMR(3,2)==1 && DMR(4,2)==0
                else
                    DMR=sortrows(DMR,1);
                end
                        
                HC=DMR(:,3);
                HS=DMR(:,4);

                clear DMR
            end
            
            %Import of data file containing DTM
            if any(volbapar==10) || any(volbapar==23)
                %If geoid undulation/height anomaly is to be computed,
                %file containg heights of the points on the Earth's surface
                %is not necessary.
                %Only DTM file with spherical harmonic coefficients is
                %necessary (computation is identical to the one in GrafLab).
            else
                loadname=get(findobj('tag','use'),'userdata');
                loadadresar=get(findobj('tag','datamat'),'userdata');

                if isempty(loadname) %Error message, if file containing heights 
                    %of the irregular surface has not been imported
                    set(findobj('tag','hlasky'),'string','Select input file with heights',...
                           'fontsize',8,'foregroundcolor','k'); drawnow;

                    warn4=warndlg('Input file containing heights of the irregular surface was not specified. Click OK and then select this file.');
                    waitfor(warn4);

                    [loadname,loadadresar]=uigetfile('*.*','Select File Containing Heights of the Irregular Surface');
                    if loadname==0
                        errordlg('Input file containing heights of the irregular surface must be specified!');

                        set(findobj('tag','hlasky'),'string','',...
                           'fontsize',8,'foregroundcolor','k'); drawnow;

                        error('Input file containing heights of the irregular surface must be specified!');
                    else
                        if find(outname=='.')>0
                            loadname=loadname(1:(find(outname=='.')-1));
                        end
                    end

                    set(findobj('tag','hlasky'),'string','',...
                           'fontsize',8,'foregroundcolor','k'); drawnow;
                end
            end

            %Entered coordinates of the grid
            fimin=str2num(get(findobj('tag','fimin'),'string')); 
            fistep=str2num(get(findobj('tag','fistep'),'string'));
            fimax=str2num(get(findobj('tag','fimax'),'string'));
            lambdamin=str2num(get(findobj('tag','lambdamin'),'string')); 
            lambdastep=str2num(get(findobj('tag','lambdastep'),'string'));
            lambdamax=str2num(get(findobj('tag','lambdamax'),'string'));
            h=str2num(get(findobj('tag','hgrid'),'string'));

            %Check of the entered coordinates
            if isempty(fimin) || isempty(fistep) || isempty(fimax) || isempty(lambdamin) || isempty(lambdastep) || isempty(lambdamax) || isempty(h) || ~isreal(fimin) || ~isreal(fistep) || ~isreal(fimax) || ~isreal(lambdamin) || ~isreal(lambdastep) || ~isreal(lambdamax) || ~isreal(h)
                errordlg('Entered grid is not correct.',...
                    'Error in point type selection');
                error('Entered grid is not correct.');       
            end

            if fimin>fimax
                errordlg('Value of Lat. min must be smaller than the Lat. max value.',...
                    'Error in irregular surface selection panel');
                error('Value of Lat. min must be smaller than the Lat. max value.'); 
            elseif fistep<=0
                errordlg('Value of Lat. step must be larger than zero.',...
                    'Error in irregular surface selection panel');
                error('Value of Lat. step must be larger than zero.'); 
            elseif lambdamin>lambdamax
                errordlg('Value of Lon. min must be smaller than Lon. max value.',...
                    'Error in irregular surface selection panel');
                error('Value of Lon. min must be smaller than Lon. max value.');
            elseif lambdastep<=0
                errordlg('Lon. step must be larger than zero.',...
                    'Error in irregular surface selection panel');
                error('Lon. step must be larger than zero.'); 
            end

            if fimin>90 || fimin<-90
                errordlg('Value of Lat. min must be within the interval <-90�,90�>.',...
                    'Error in irregular surface selection panel');
                error('Value of Lat. min must be within the interval <-90�,90�>.');
            end
            if fimax>90 || fimax<-90
                errordlg('Value of Lat. max must be within the interval <-90�,90�>.',...
                    'Error in irregular surface selection panel');
                error('Value of Lat. max must be within the interval <-90�,90�>.');
            end
            if lambdamin>360 || lambdamin<-180
                errordlg('Value of Lon. min must be within the interval <-180�,180�> or <0�,360�>.',...
                    'Error in irregular surface selection panel');
                error('Value of Lon. min must be within the interval <-180�,180�> or <0�,360�>.');
            end
            if lambdamax>360 || lambdamax<-180
                errordlg('Value of Lon. max must be within the interval <-180�,180�> or <0�,360�>.',...
                    'Error in irregular surface selection panel');
                error('Value of Lon. max must be within the interval <-180�,180�> or <0�,360�>.');
            end
            if (lambdamax-lambdamin)>360
                errordlg('Longitude must be in the range <-180�,180�> or <0�,360�>.',...
                    'Error in irregular surface selection panel');
                error('Longitude must be in the range <-180�,180�> or <0�,360�>.');
            end
            
            %Vectors phi and lambda in one longitude and latitude
            %parallel, respectively
            fi=(fimin:fistep:(fimax+1000*eps))';
            lambda=(lambdamin:lambdastep:(lambdamax+1000*eps))';
            %The value "1000*eps" added the to maximum values "fimax" and
            %"lambdamax" is used in order to avoid a possible rounding erros 
            %in defining the grid, especially if the defining coordinates
            %of the grid have a lot of digits. Without the value "1000*eps" 
            %the number of points in the defined grid would not have to be 
            %the same (because of the rounding errors) as the number of 
            %points in the file containing heights of the irregular surface.
            %The value "x = 1000*eps" can be changed. However, one should 
            %keep in mind that the inequalities "fistep>x" and "lambdastep>x" 
            %must hold true.
            fi=pi/180*(fi(:));
            lambda=pi/180*(lambda(:));

            %Grid, which is to be displayed has to have at least two
            %points in latitude parallels and at least two points in
            %longitude parallels.
            if display_data==1
                if length(fi)<2 || length(lambda)<2
                    warn3=warndlg('To display computed data on a grid, it must contain at least two distinct points in one parallel and two distinct points in one meridian. After clicking OK, the computation will start, although the data will not be displayed.');
                    waitfor(warn3);
                    display_data=0;
                end
            end
            
            if coord==1 %Entered spherical coordinates                       
                
                %Spherical radius
                r=(R+h)*ones(length(fi),1);
                hsph=h;
                
                %Spherical latitude
                fiG=fi;  
                
                %Transformation of spherical latitude into the ellipsoidal
                %latitude
                [X,Y,Z]=sph2cart(0*fiG,fiG,r);
                [fi, lambda_del, h]=ecef2geodetic(X,Y,Z,[aEl eEl]);

                clear X Y Z lambda_del
            elseif coord==0 %Entered ellipsoidal coordinates
                
                %Trasformation of (fi, lambda, h) into (X, Y, Z)
                [X,Y,Z]=geodetic2ecef(fi,0*zeros(length(fi),1),h*ones(length(fi),1),[aEl eEl]);  
                r=sqrt(X.*X+Y.*Y+Z.*Z); %Radius

                %Spherical latitude
                fiG=atan(Z./sqrt(X.*X+Y.*Y)); 
                
                clear X Y Z 
            end

            %Computation of the coefficients C0,0; C2,0; ...; C20,0
            %of the selected ellipsoid
            CEl=zeros(length(C),1);
            for n=0:10
                CEl(2*n==stupen & rad==0,1)=((-1)^n*(3*eEl^(2*n))/((2*n+1)*(2*n+3)*sqrt(4*n+1))*(1-n-5^(3/2)*n*CEl_20/eEl^2)).*(aEl./R).^(2*n).*(GMEl/GM);
            end                                
     
            if any(volbapar==11) || any(volbapar==12) || any(volbapar==13) || any(volbapar==14) || any(volbapar==15) || any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20) || any(volbapar==25)
                grav=1;
            else
                grav=0;
            end

            if  any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==5) || any(volbapar==6) || any(volbapar==7) || any(volbapar==8) || any(volbapar==9) || any(volbapar==10) || any(volbapar==19) || any(volbapar==21) || any(volbapar==22) || any(volbapar==23) || any(volbapar==24)
                por=1;
                deltaC=C-CEl;
            else
                por=0;
            end
 
            if any(volbapar==20)
                normal=1;
            else
                normal=0;
            end

            clear GGM stupen rad                   
            if normal==0
                clear CEl
            end
            
            %Initialization for geoid undulation/height anomaly
            %computation
            N1c=0; H=0;

            %If geoid/height anomaly is to be computed
            if any(volbapar==10) || any(volbapar==23)
                geoid=1;
                if h~=0
                    errordlg('To compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                        'Error in irregular surface selection panel');
                    error('To compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                end
            else
                geoid=0;
            end

            %Indices of the spherical harmonic coefficients
            index=zeros(nmax+1,1);
            index(1,1)=1;
            for i=1:nmax
                index(i+1,1)=index(i,1)+i;
            end

            %Initialization of the matrices and vectors for the computation of fnALFs
            Pnm=zeros(length(fiG),nmax+1);
            q=(R./r);
            q2=(R./r).^2;
            u=cos(fiG);
            t=sin(fiG);
            
            %Initialization for extended-range arithmetic approach
            if volbaALFs==3
                                       
                bit=mexext; %Bit version of Matlab
                bit=bit(end-1:end);
                bit=str2double(bit);
                if bit==32
                    bit=32;
                elseif bit==64
                    bit=64;
                else
                    bit=64;
                end
                
                nmax23=nmax*2+3;
                rr=zeros(nmax23,1); ri=rr;
                dd=zeros(nmax,1); am=dd; bm=am;

                m1=1:nmax23;
                rr(m1)=sqrt(m1);
                ri(m1)=1./rr;
                m2=1:nmax;
                dd(m2)=rr(2*m2+3).*ri(2*m2+2);

                IND=960;
                BIG=2^IND;
                BIGI=2^(-IND);
                BIGS=2^(IND/2);
                BIGSI=2^(-IND/2);
                ROOT3=1.732050807568877;
                
                if bit==32
                    pm=am;
                    ps1=zeros(length(fiG),nmax); 
                    ips1=ps1;
                    x=ROOT3*u.*q;
                    ix=zeros(size(x));
                    ps1(:,1)=x;
                    ips1(:,1)=ix;
                    for m3=2:nmax
                        x=(dd(m3-1)*u).*x.*q;
                        y=abs(x);
                        iy=y>=BIGS;
                        if any(iy)
                            x(iy)=x(iy)*BIGI;
                            ix(iy)=ix(iy)+1;
                        end
                        iy=y<BIGSI;
                        if any(iy)
                            x(iy)=x(iy)*BIG;
                            ix(iy)=ix(iy)-1;
                        end
                        ps1(:,m3)=x;
                        ips1(:,m3)=ix;
                    end
                elseif bit==64
					tq=t.*q;
                    temp1=zeros(length(fiG),1);
                    temp2=ones(length(fiG),1);
                    temp3=temp2;
                    temp4=temp1;
                    temp5=temp1+BIGI;
                    ps1b=zeros(length(fiG),nmax); 
                    ips1b=ps1b;
                    xb=ROOT3*u.*q;
                    ixb=zeros(size(xb));
                    ps1b(:,1)=xb;
                    ips1b(:,1)=ixb;
                    for m3=2:nmax
                        xb=(dd(m3-1)*u).*xb.*q;
                        yb=abs(xb);
                        iyb=yb>=BIGS;
                        if any(iyb)
                            xb(iyb)=xb(iyb)*BIGI;
                            ixb(iyb)=ixb(iyb)+1;
                        end
                        iyb=yb<BIGSI;
                        if any(iyb)
                            xb(iyb)=xb(iyb)*BIG;
                            ixb(iyb)=ixb(iyb)-1;
                        end
                        ps1b(:,m3)=xb;
                        ips1b(:,m3)=ixb;
                    end
                end
                    
                clear dd
            end

            %Initialization of the matrices and vectors for the 
            %computation of the first-order derivatives of fnALFs
            if any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16)  || any(volbapar==20)
                dALFs=1;
                dPnm=zeros(length(fi),nmax+1);
                qu=q./u;
                tu=t./u;
                
                %Treatment of the dPnm singularity
                singdPnm=fi==pi/2 | fi==-pi/2;
            else
                dALFs=0;
            end   
            
            %Initialization of the matrices and vectors for the 
            %computation of the second-order derivatives of fnALFs
            if any(volbapar==6) || any(volbapar==12)
                ddALFs=1;
                ddPnm=zeros(length(fi),nmax+1);
                
                %Treatment of the ddPnm singularity
                singddPnm=fi==pi/2 | fi==-pi/2;
            else
                ddALFs=0;
            end   
            
            %Status line
            progressbar=findobj('tag','hlasky');

            %% Summation over m
            for m=nmax:-1:0

                %Update of the progress bar
                if rem(m,10)==0
                    set(progressbar,'string',...
                        sprintf('Progress: m = %5.0d',m),...
                        'fontsize',8); drawnow;
                end

                %Selection of the spherical harmonic coefficients of order m
                %======================================================
                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                    Cm=C(index((m+1):end)+m);
                end

                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                    deltaCm=deltaC(index((m+1):end)+m);
                end
                                
                if normal==1 %C's spherical harmonic coefficients for the functionals with the normal field
                    if m==0
                        CElm=CEl(index((m+1):end)+m);
                    end
                end
                            
                if geoid==1
                    HCm=HC(index((m+1):end)+m); 
                    HSm=HS(index((m+1):end)+m);
                end
                
                Sm=S(index((m+1):end)+m);
                %======================================================

                %% Computation of the modified fnALFs
                if volbaALFs==1 %Standard forward column method
                    if m==0
                        Pnm(:,1)=1;
                    elseif m==1                    
                        Pnm(:,1)=sqrt(3)*u.*q;  
                    elseif m>1                            
                        i=2*(2:m);
                        i1=sqrt((i+ones(size(i)))./i);
                        Pnm(:,1)=u.^m*sqrt(3)*prod(i1).*q.^m;
                    end

                    if m==nmax
                    elseif m<=(nmax-1)
                        n=m+1;
                        anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                        Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                    end

                    if m<(nmax-1)
                        j=3;
                        for n=m+2:nmax                            
                            anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                            bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                            Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                            j=j+1;
                        end
                    end

                elseif volbaALFs==2 %Modified forward column method
                    if m==0
                        Pnm(:,1)=1e-280;
                    elseif m==1

                        Pnm(:,1)=sqrt(3)*q*1e-280;  
                    elseif m>1                            
                        i=2*(2:m);
                        i1=sqrt((i+ones(size(i)))./i);
                        Pnm(:,1)=sqrt(3)*prod(i1)*(q.^m)*1e-280;
                    end

                    if m==nmax
                    elseif m<=(nmax-1)
                        n=m+1;
                        anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                        Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                    end

                    if m<(nmax-1)
                        j=3;
                        for n=m+2:nmax                            
                            anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                            bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                            Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                            j=j+1;
                        end
                    end

                elseif volbaALFs==3 %Extended-range arithmetic 
                    if bit==32 %32 bit version of Matlab
             
                        am(m+1)=rr(2*m+3);
                        for n=m+2:nmax
                            ww=rr(2*n+1)*ri(n-m)*ri(n+m);
                            am(n)=rr(2*n-1)*ww;
                            bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*ww;
                        end

                        if m~=0
                            for i=1:length(fiG) 
                                x=ps1(i,m);
                                ix=ips1(i,m);

                                if ix==0
                                    pm(m)=x;
                                elseif ix<-1
                                    pm(m)=0;  
                                elseif ix<0
                                    pm(m)=x*BIGI;
                                else
                                    pm(m)=x*BIG;
                                end

                                if m==nmax
                                    Pnm(i,1:(nmax-m+1))=pm(m:end);
                                    continue;
                                end

                                y=x;
                                iy=ix;
                                x=(am(m+1)*t(i)*q(i))*y;
                                ix=iy;
                                w=abs(x);

                                if w>=BIGS
                                    x=x*BIGI;
                                    ix=ix+1;
                                elseif w<BIGSI
                                    x=x*BIG;
                                    ix=ix-1;
                                end

                                if ix==0
                                    pm(m+1)=x;
                                elseif ix<-1  
                                    pm(m+1)=0;    
                                elseif ix<0
                                    pm(m+1)=x*BIGI;
                                else
                                    pm(m+1)=x*BIG;
                                end

                                for n=m+2:nmax 
                                    id=ix-iy;

                                    if id==0
                                        zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*y;
                                        iz=ix;
                                    elseif id==1
                                        zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*(y*BIGI);
                                        iz=ix;
                                    elseif id==-1
                                        zz=(am(n)*t(i)*q(i))*(x*BIGI)-bm(n)*q2(i)*y;
                                        iz=iy;
                                    elseif id>1
                                        zz=(am(n)*t(i)*q(i))*x;
                                        iz=ix;
                                    else
                                        zz=-bm(n)*q2(i)*y;
                                        iz=iy;
                                    end

                                    w=abs(zz);

                                    if w>=BIGS
                                        zz=zz*BIGI;
                                        iz=iz+1;
                                    elseif w<BIGSI
                                        zz=zz*BIG;
                                        iz=iz-1;
                                    end

                                    if iz==0
                                        pm(n)=zz;
                                    elseif iz<-1
                                        pm(n)=0;     
                                    elseif iz<0
                                        pm(n)=zz*BIGI;
                                    else
                                        pm(n)=zz*BIG;
                                    end

                                    y=x;
                                    iy=ix;
                                    x=zz;
                                    ix=iz;                                           
                                end

                                Pnm(i,1:(nmax-m+1))=pm(m:end);  
                            end

                        elseif m==0
                            Pnm(:,1)=1;
                            Pnm(:,2)=sqrt(3)*t.*q;

                            for i=2:nmax
                                Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                            end

                            clear rr ri am bm pm ps1 ips1 m1 m2 ...
                                dd ix x y iy w iz zz
                        end
                    elseif bit==64 %64 bit version of Matlab
                        
                        am(m+1)=rr(2*m+3);
                        for n=m+2:nmax
                            ww=rr(2*n+1)*ri(n-m)*ri(n+m);
                            am(n)=rr(2*n-1)*ww;
                            bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*ww;
                        end

                        if m==0 %Zonal modified fnALFs
                            Pnm(:,1)=1;
                            Pnm(:,2)=sqrt(3)*tq;

                            for i=2:nmax
                                Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*tq-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                            end

                            clear rr ri am bm ps1b ips1b m1 m2 dd ...
                                ixb xb yb iyb wb izb zzb pmxb pm0b pmxBIGIb ...
                                pmxBIGb wb wBIGSb wBIGSIb pm1xb pm10b ...
                                pm1xBIGIb pm1xBIGb idb id0b id1b id_1b ...
                                idv1b idm1b iz0b izm_1b izm0b izv0b tq temp1...
                                temp2 temp3 temp4 temp5

                        elseif m~=0 %Non-zonal modified fnALFs                                    
                            xb=ps1b(:,m);
                            ixb=ips1b(:,m);
                                                               
                            temp5(ixb==0)=1;
                            temp5(ixb<-1)=0;
                            %temp5(izb>=-1 & izb<0)=BIGI;
                            %The condition "izb>=-1 & izb<0"
                            %is useless, as "izb" is already
                            %initialized as "izb=BIGI".
                            temp5(ixb>0)=BIG;
                                   
                            Pnm(:,1)=xb.*temp5;
                            temp5=temp5.*0+BIGI; 
                                                              
                            if m<nmax
                               yb=xb;
                               iyb=ixb;

                               xb=(am(m+1).*tq).*yb;
                               ixb=iyb;
                               wb=abs(xb);

                               wBIGSb=wb>=BIGS;
                               wBIGSIb=wb<BIGSI;
                               temp3(wBIGSb)=BIGI;
                               temp3(wBIGSIb)=BIG;
                               temp4(wBIGSb)=1;
                               temp4(wBIGSIb)=-1;

                               xb=xb.*temp3;
                               ixb=ixb+temp4;
                               temp3=temp2;
                               temp4=temp4.*0;
                               
                               temp5(ixb==0)=1;
                               temp5(ixb<-1)=0;
                               %temp5(izb>=-1 & izb<0)=BIGI;
                               %The condition "izb>=-1 & izb<0"
                               %is useless, as "izb" is already
                               %initialized as "izb=BIGI".
                               temp5(ixb>0)=BIG;
                                   
                               Pnm(:,2)=xb.*temp5;
                               temp5=temp5.*0+BIGI; 

                               for n=m+2:nmax
                                   idb=ixb-iyb;

                                   id0b=idb==0;
                                   id1b=idb==1;
                                   id_1b=idb==-1;
                                   idv1b=idb>1;
                                   
                                   temp1(id0b)=1;
                                   temp1(id1b)=1;
                                   temp2(id1b)=BIGI;
                                   temp1(id_1b)=BIGI;
                                   temp1(idv1b)=1;
                                   temp2(idv1b)=0;
                                   
                                   zzb=(am(n).*tq).*(xb.*temp1)-bm(n).*((yb.*q2).*temp2);
                                   izb=iyb;
                                   id0b_id1b_idv1b=id0b | id1b | idv1b;
                                   izb(id0b_id1b_idv1b)=ixb(id0b_id1b_idv1b);
                                   temp1=temp1.*0;
                                   temp2=temp1+1;

                                   wb=abs(zzb);

                                   wBIGSb=wb>=BIGS;
                                   wBIGSIb=wb<BIGSI;
                                   temp3(wBIGSb)=BIGI;
                                   temp3(wBIGSIb)=BIG;
                                   temp4(wBIGSb)=1;
                                   temp4(wBIGSIb)=-1;

                                   zzb=zzb.*temp3;
                                   izb=izb+temp4;
                                   temp3=temp2;
                                   temp4=temp4.*0;

                                   temp5(izb==0)=1;
                                   temp5(izb<-1)=0;
                                   %temp5(izb>=-1 & izb<0)=BIGI;
                                   %The condition "izb>=-1 & izb<0"
                                   %is useless, as "izb" is already
                                   %initialized as "izb=BIGI".
                                   temp5(izb>0)=BIG;
                                   
                                   Pnm(:,n-m+1)=zzb.*temp5;
                                   temp5=temp1+BIGI;   
                        
                                   yb=xb;
                                   iyb=ixb;
                                   xb=zzb;
                                   ixb=izb;
                               end   
                            end
                        end	
                    end
                end
                
                %If nmin~=0
                %======================================================
                if nmin~=0
                    if LNOFnmin==0
                        if por==1
                            deltaCm(1:(nmin-m))=0;
                        end

                        if grav==1
                            Cm(1:(nmin-m))=0;
                        end

                        if geoid==1
                            HCm(1:(nmin-m))=0;
                            HSm(1:(nmin-m))=0;
                        end

                        Sm(1:(nmin-m))=0;
                    elseif LNOFnmin==1
                        Pnm(:,1:(nmin-m))=0;
                    end
                end
                %======================================================
   
                %% Computation of the first-order derivatives of the modified fnALFs
                if dALFs==1  
                    if volbaALFs==1 || volbaALFs==3
                        enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                        if m==0 %Zonal modified dALFs
                            dPnm(:,1)=0.*u;
                            dPnm(:,2)=sqrt(3)*u.*q;
                            dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                        elseif m==nmax %Sectorial modified dALFs
                            dPnm(:,1)=-m*(tu).*Pnm(:,1);
                        else
                            dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                            dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne ALFs
                        end

                    elseif volbaALFs==2
                        enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                        if m==0 %Zonal modified dALFs
                            dPnm(:,1)=0.*u*1e-280;
                            dPnm(:,2)=sqrt(3)*u.*q*1e-280;
                            dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                        elseif m==nmax
                            dPnm(:,1)=-m*(tu).*Pnm(:,1); %Sectorial modified dALFs
                        else
                            dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                            dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne dALFs
                        end
                    end

                    %Treatment of the dALFs singularity
                    dPnm(singdPnm,:)=0;

                    if ddALFs==1 %If the second-order derivatives of the modified fnALFs are to be computed
                        
                        if m==0 %Zonal modified ddALFs
                            ddPnm=bsxfun(@times,tu,dPnm)-bsxfun(@times,(0:nmax).*((0:nmax)+1),Pnm);
                        else
                            ddPnm(:,1:end-m)=bsxfun(@times,tu,dPnm(:,1:end-m))+bsxfun(@times,m^2./u.^2,Pnm(:,1:end-m))-bsxfun(@times,(m:nmax).*((m:nmax)+1),Pnm(:,1:end-m));
                        end
                                                        
                        %Treatment of the ddALFs singularity
                        ddPnm(singddPnm,:)=0;
                    end                                                                               
                end

                %% Loop for 1:NF (number of computing functionals)
                for i=1:pocetpar                   
                    if volbapar(i)==1       
                    elseif volbapar(i)==2 %Deflection of the vertical eta                       
                        
                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_eta{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_eta{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_eta{K+1}=[ampl_eta{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end 
                                end
                                
                                ampl_eta{K+1}=strrep(ampl_eta{K+1},')(',').*(');
                                ampl_eta{K+1}=eval(ampl_eta{K+1});
                            end
                            
                            Lm=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_eta{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3 %Computation using the standard 
                                %forward column method or extended-range arithmetic
                                if m==nmax
                                    Aeta{K+1}=zeros(length(fiG),nmax+1);
                                    Beta{K+1}=zeros(length(fiG),nmax+1);
                                end

                                Aeta{K+1}(:,m+1)=Lm*deltaCm;
                                Beta{K+1}(:,m+1)=Lm*Sm;                        
                            elseif volbaALFs==2 %Computation using the modified forward column method combined with Horner's scheme
                                if m==nmax
                                    eta{K+1}=zeros(length(fiG),length(lambda));
                                end

                                eta{K+1}=bsxfun(@times,eta{K+1},u)+(-Lm*deltaCm*sin(m*lambda')+Lm*Sm*cos(m*lambda'));
                            end
                        end
                        
                    elseif volbapar(i)==3 %Deflection of the vertical xi

                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_ksi{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_ksi{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_ksi{K+1}=[ampl_ksi{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end                                        
                                end
                                
                                ampl_ksi{K+1}=strrep(ampl_ksi{K+1},')(',').*(');
                                ampl_ksi{K+1}=eval(ampl_ksi{K+1});
                            end
                            
                            dLm=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_ksi{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    Aksi{K+1}=zeros(length(fiG),nmax+1);
                                    Bksi{K+1}=zeros(length(fiG),nmax+1);
                                end

                                Aksi{K+1}(:,m+1)=dLm*deltaCm;
                                Bksi{K+1}(:,m+1)=dLm*Sm; 
                            elseif volbaALFs==2
                                if m==nmax
                                    ksi{K+1}=zeros(length(fiG),length(lambda));
                                end

                                ksi{K+1}=bsxfun(@times,ksi{K+1},u)+(dLm*deltaCm*cos(m*lambda')+dLm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==4 %Deflection of the vertical Theta

                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_Tksi{K+1}='ones(1,nmax+1)';
                                    ampl_Teta{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Tksi{K+1}='ones(1,nmax+1)';
                                    ampl_Teta{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Tksi{K+1}=[ampl_Tksi{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                        ampl_Teta{K+1}=[ampl_Teta{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end  
                                end
                                
                                ampl_Tksi{K+1}=strrep(ampl_Tksi{K+1},')(',').*(');
                                ampl_Tksi{K+1}=eval(ampl_Tksi{K+1});
                                ampl_Teta{K+1}=strrep(ampl_Teta{K+1},')(',').*(');
                                ampl_Teta{K+1}=eval(ampl_Teta{K+1});
                            end

                            Lm=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_Teta{K+1}(:,(m+1):end));
                            dLm=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_Tksi{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    ATeta{K+1}=zeros(length(fiG),nmax+1);
                                    BTeta{K+1}=zeros(length(fiG),nmax+1);
                                    ATksi{K+1}=ATeta{K+1};
                                    BTksi{K+1}=BTeta{K+1};
                                end

                                ATeta{K+1}(:,m+1)=Lm*deltaCm;
                                BTeta{K+1}(:,m+1)=Lm*Sm;
                                ATksi{K+1}(:,m+1)=dLm*deltaCm;
                                BTksi{K+1}(:,m+1)=dLm*Sm; 
                            elseif volbaALFs==2
                                if m==nmax
                                    Teta{K+1}=zeros(length(fiG),length(lambda));
                                    Tksi{K+1}=zeros(length(fiG),length(lambda));
                                end

                                Teta{K+1}=bsxfun(@times,Teta{K+1},u)+(-Lm*deltaCm*sin(m*lambda')+Lm*Sm*cos(m*lambda'));
                                Tksi{K+1}=bsxfun(@times,Tksi{K+1},u)+(dLm*deltaCm*cos(m*lambda')+dLm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==5 %Disturbing potential

                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_T{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_T{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_T{K+1}=[ampl_T{K+1} sprintf('((0:nmax)+%d)',KK+1)];
                                    end
                                end
                                
                                ampl_T{K+1}=strrep(ampl_T{K+1},')(',').*(');
                                ampl_T{K+1}=eval(ampl_T{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_T{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AT{K+1}=zeros(length(fiG),nmax+1);
                                    BT{K+1}=zeros(length(fiG),nmax+1);
                                end

                                AT{K+1}(:,m+1)=Lm*deltaCm;
                                BT{K+1}(:,m+1)=Lm*Sm;    
                            elseif volbaALFs==2
                                if m==nmax
                                    T{K+1}=zeros(length(fiG),length(lambda));
                                end

                                T{K+1}=bsxfun(@times,T{K+1},u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end
                        
                    elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                        
                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Trr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Tppll{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Trr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Tppll{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Trr{K+1}=[ampl_Trr{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                        ampl_Tppll{K+1}=[ampl_Tppll{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Trr{K+1}=strrep(ampl_Trr{K+1},')(',').*(');
                                ampl_Tppll{K+1}=strrep(ampl_Tppll{K+1},')(',').*(');
                                ampl_Trr{K+1}=eval(ampl_Trr{K+1});
                                ampl_Tppll{K+1}=eval(ampl_Tppll{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Trr{K+1}(:,(m+1):end));
                            Lmff=bsxfun(@times,ddPnm(:,1:(nmax-m+1)),ampl_Tppll{K+1}(:,(m+1):end));
                            Lmll=bsxfun(@times,m^2*Pnm(:,1:(nmax-m+1)),ampl_Tppll{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    ATrr{K+1}=zeros(length(fiG),nmax+1);
                                    BTrr{K+1}=zeros(length(fiG),nmax+1);
                                    ATff{K+1}=ATrr{K+1};
                                    BTff{K+1}=ATrr{K+1};
                                    ATll{K+1}=ATrr{K+1};
                                    BTll{K+1}=ATrr{K+1};
                                end

                                ATrr{K+1}(:,m+1)=Lm*deltaCm;
                                BTrr{K+1}(:,m+1)=Lm*Sm;  

                                ATff{K+1}(:,m+1)=Lmff*deltaCm;
                                BTff{K+1}(:,m+1)=Lmff*Sm;

                                ATll{K+1}(:,m+1)=Lmll*deltaCm;
                                BTll{K+1}(:,m+1)=Lmll*Sm;
                            elseif volbaALFs==2
                                if m==nmax
                                    Trr{K+1}=zeros(length(fiG),length(lambda));
                                    Tff{K+1}=Trr{K+1};
                                    Tll{K+1}=Trr{K+1};
                                end

                                Trr{K+1}=bsxfun(@times,Trr{K+1},u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                Tff{K+1}=bsxfun(@times,Tff{K+1},u)+(Lmff*deltaCm*cos(m*lambda')+Lmff*Sm*sin(m*lambda'));
                                Tll{K+1}=bsxfun(@times,Tll{K+1},u)+(Lmll*deltaCm*cos(m*lambda')+Lmll*Sm*sin(m*lambda'));
                            end
                        end
                        
                    elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                        
                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Trf{K+1}='((0:nmax)+1)';
                                    ampl_Tfl{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Trf{K+1}='((0:nmax)+1)';
                                    ampl_Tfl{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Trf{K+1}=[ampl_Trf{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                        ampl_Tfl{K+1}=[ampl_Tfl{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Trf{K+1}=strrep(ampl_Trf{K+1},')(',').*(');
                                ampl_Tfl{K+1}=strrep(ampl_Tfl{K+1},')(',').*(');
                                ampl_Trf{K+1}=eval(ampl_Trf{K+1});
                                ampl_Tfl{K+1}=eval(ampl_Tfl{K+1});
                            end

                            Lmrf=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_Trf{K+1}(:,(m+1):end));
                            Lmrl=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_Trf{K+1}(:,(m+1):end));
                            Lmfl=bsxfun(@times,m*dPnm(:,1:(nmax-m+1)),ampl_Tfl{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    ATrf{K+1}=zeros(length(fiG),nmax+1);
                                    BTrf{K+1}=zeros(length(fiG),nmax+1);
                                    ATrl{K+1}=ATrf{K+1};
                                    BTrl{K+1}=ATrf{K+1};
                                    ATfl{K+1}=ATrf{K+1};
                                    BTfl{K+1}=ATrf{K+1};
                                end

                                ATrf{K+1}(:,m+1)=Lmrf*deltaCm;
                                BTrf{K+1}(:,m+1)=Lmrf*Sm;  

                                ATrl{K+1}(:,m+1)=Lmrl*deltaCm;
                                BTrl{K+1}(:,m+1)=Lmrl*Sm;

                                ATfl{K+1}(:,m+1)=Lmfl*deltaCm;
                                BTfl{K+1}(:,m+1)=Lmfl*Sm;
                            elseif volbaALFs==2
                                if m==nmax
                                    Trf{K+1}=zeros(length(fiG),length(lambda));
                                    Trl{K+1}=Trf{K+1};
                                    Tfl{K+1}=Trf{K+1};
                                end

                                Trf{K+1}=bsxfun(@times,Trf{K+1},u)+(Lmrf*deltaCm*cos(m*lambda')+Lmrf*Sm*sin(m*lambda'));
                                Trl{K+1}=bsxfun(@times,Trl{K+1},u)+(Lmrl*Sm*cos(m*lambda')-Lmrl*deltaCm*sin(m*lambda'));
                                Tfl{K+1}=bsxfun(@times,Tfl{K+1},u)+(Lmfl*Sm*cos(m*lambda')-Lmfl*deltaCm*sin(m*lambda'));
                            end
                        end
                        
                    elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                        
                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Tzz{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Txxyy{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Tzz{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Txxyy{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Tzz{K+1}=[ampl_Tzz{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                        ampl_Txxyy{K+1}=[ampl_Txxyy{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Tzz{K+1}=strrep(ampl_Tzz{K+1},')(',').*(');
                                ampl_Txxyy{K+1}=strrep(ampl_Txxyy{K+1},')(',').*(');
                                ampl_Tzz{K+1}=eval(ampl_Tzz{K+1});
                                ampl_Txxyy{K+1}=eval(ampl_Txxyy{K+1});
                            end

                            if m<2
                                Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Tzz{K+1}(:,(m+1):end));

                                %Txx
                                bnm=((m:nmax)+m+1).*((m:nmax)+m+2)./2./(m+1);
                                Pnmxxyy=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Txxyy{K+1}(:,(m+1):end));
                                LmTxx1=bsxfun(@times,Pnmxxyy,bnm-((m:nmax)+1).*((m:nmax)+2));
                                LmTyy1=bsxfun(@times,Pnmxxyy,bnm);

                                if m==0
                                    anm=sqrt(2)/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                else
                                    anm=1./4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                end

                                LmTxx2=bsxfun(@times,Pnmxxyy,anm);
                                LmTxx3=zeros(length(fi),nmax-m+1); %Coefficients cnm are equal to zeros, if m==0 a m==1
                            else
                                Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Tzz{K+1}(:,(m+1):end));

                                %Txx
                                bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                Pnmxxyy=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Txxyy{K+1}(:,(m+1):end));
                                LmTxx1=bsxfun(@times,Pnmxxyy,bnm-((m:nmax)+1).*((m:nmax)+2));
                                LmTyy1=bsxfun(@times,Pnmxxyy,bnm);

                                anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                LmTxx2=bsxfun(@times,Pnmxxyy,anm);

                                if m==2
                                    cnm=sqrt(2)/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                else
                                    cnm=1./4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                end

                                LmTxx3=bsxfun(@times,Pnmxxyy,cnm);                                    
                            end

                            if m==nmax
                                ATzz{K+1}=zeros(length(fiG),nmax+1);
                                BTzz{K+1}=zeros(length(fiG),nmax+1);
                                ATxx1{K+1}=ATzz{K+1};
                                BTxx1{K+1}=ATzz{K+1};
                                ATxx2{K+1}=ATzz{K+1};
                                BTxx2{K+1}=ATzz{K+1};
                                ATxx3{K+1}=ATzz{K+1};
                                BTxx3{K+1}=BTzz{K+1};
                                ATyy1{K+1}=ATzz{K+1};
                                BTyy1{K+1}=BTzz{K+1};
                            end

                            ATzz{K+1}(:,m+1)=Lm*deltaCm;
                            BTzz{K+1}(:,m+1)=Lm*Sm;  

                            %Txx
                            ATxx1{K+1}(:,m+1)=LmTxx1*deltaCm;
                            BTxx1{K+1}(:,m+1)=LmTxx1*Sm;

                            if m<2
                                ATxx2{K+1}(:,m+1)=LmTxx2*deltaC(index((m+1):end)+m+2);
                                BTxx2{K+1}(:,m+1)=LmTxx2*S(index((m+1):end)+m+2);

                                ATxx3{K+1}(:,m+1)=0;
                                BTxx3{K+1}(:,m+1)=0;
                            elseif m>nmax-2
                                ATxx2{K+1}(:,m+1)=0;
                                BTxx2{K+1}(:,m+1)=0;

                                ATxx3{K+1}(:,m+1)=LmTxx3*deltaC(index((m+1):end)+m-2);
                                BTxx3{K+1}(:,m+1)=LmTxx3*S(index((m+1):end)+m-2);
                            else
                                ATxx2{K+1}(:,m+1)=LmTxx2*deltaC(index((m+1):end)+m+2);
                                BTxx2{K+1}(:,m+1)=LmTxx2*S(index((m+1):end)+m+2);

                                ATxx3{K+1}(:,m+1)=LmTxx3*deltaC(index((m+1):end)+m-2);
                                BTxx3{K+1}(:,m+1)=LmTxx3*S(index((m+1):end)+m-2);
                            end

                            %Tyy
                            ATyy1{K+1}(:,m+1)=LmTyy1*deltaCm;
                            BTyy1{K+1}(:,m+1)=LmTyy1*Sm;
                        end
                        
                        %Modified forward column method cannot be
                        %applied 
                        
                    elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                               
                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Tnondiag{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Tnondiag{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Tnondiag{K+1}=[ampl_Tnondiag{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Tnondiag{K+1}=strrep(ampl_Tnondiag{K+1},')(',').*(');
                                ampl_Tnondiag{K+1}=eval(ampl_Tnondiag{K+1});
                            end
                            
                            if m<2                                       
                                if m==0
                                    betanm=((m:nmax)+2)./2.*sqrt(1+ones(1,nmax+1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=zeros(1,nmax+1);

                                    %Txy
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);  
                                    
                                    Pnmnondiag=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ampl_Tnondiag{K+1}(:,(m+1):end));
                                    LmTxy1=bsxfun(@times,Pnmnondiag,dnm);
                                    LmTxy2=zeros(length(fi),nmax+1);
                                    LmTxy3=zeros(length(fi),nmax+1);

                                    %Tyz
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    LmTyz1=bsxfun(@times,Pnmnondiag,minm);
                                    LmTyz2=zeros(length(fi),nmax+1);
                                else
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2).*sqrt((m:nmax).*((m:nmax)+1)./2);

                                    %Txy
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);  
                                    
                                    Pnmnondiag=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ampl_Tnondiag{K+1}(:,(m+1):end));
                                    LmTxy1=bsxfun(@times,Pnmnondiag,dnm);

                                    gnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+1).*sqrt((m:nmax)-1).*((m:nmax)+2);
                                    LmTxy2=bsxfun(@times,Pnmnondiag,gnm);

                                    LmTxy3=zeros(length(fi),nmax);

                                    %Tyz
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    LmTyz1=bsxfun(@times,Pnmnondiag,minm);

                                    LmTyz2=zeros(length(fi),nmax);
                                end

                                %Txz
                                Pnmnondiag=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Tnondiag{K+1}(:,(m+1):end));
                                LmTxz1=bsxfun(@times,Pnmnondiag,betanm);
                                LmTxz2=bsxfun(@times,Pnmnondiag,gamanm); 
                            else
                                Pnmnondiag=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Tnondiag{K+1}(:,(m+1):end));
                                
                                %Txz
                                betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);

                                LmTxz1=bsxfun(@times,Pnmnondiag,betanm);
                                LmTxz2=bsxfun(@times,Pnmnondiag,gamanm);
                                
                                Pnmnondiag=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ampl_Tnondiag{K+1}(:,(m+1):end));

                                %Txy
                                dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                LmTxy1=bsxfun(@times,Pnmnondiag,dnm);

                                gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                LmTxy2=bsxfun(@times,Pnmnondiag,gnm);

                                if m==2
                                    %Txy
                                    LmTxy3=zeros(length(fi),nmax-1);    
                                elseif m==3
                                    %Txy
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                    LmTxy3=bsxfun(@times,Pnmnondiag,hnm);
                                else
                                    %Txy
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                    LmTxy3=bsxfun(@times,Pnmnondiag,hnm);
                                end

                                %Tyz
                                minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                LmTyz1=bsxfun(@times,Pnmnondiag,minm);

                                ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                LmTyz2=bsxfun(@times,Pnmnondiag,ninm);
                            end

                            if m==nmax
                                ATxz1{K+1}=zeros(length(fiG),nmax+1);
                                BTxz1{K+1}=zeros(length(fiG),nmax+1);
                                ATxz2{K+1}=ATxz1{K+1};
                                BTxz2{K+1}=ATxz1{K+1};
                                ATxy1{K+1}=ATxz1{K+1};
                                BTxy1{K+1}=ATxz1{K+1};
                                ATxy2{K+1}=ATxz1{K+1};
                                BTxy2{K+1}=ATxz1{K+1};
                                ATxy3{K+1}=ATxz1{K+1};
                                BTxy3{K+1}=ATxz1{K+1};
                                ATyz1{K+1}=ATxz1{K+1};
                                BTyz1{K+1}=ATxz1{K+1};
                                ATyz2{K+1}=ATxz1{K+1};
                                BTyz2{K+1}=ATxz1{K+1};
                            end

                            ATxy2{K+1}(:,m+1)=(LmTxy2*deltaCm).*q;
                            BTxy2{K+1}(:,m+1)=(LmTxy2*Sm).*q;

                            if m<2
                                ATxz1{K+1}(:,m+1)=LmTxz1*deltaC(index((m+1):end)+m+1);
                                BTxz1{K+1}(:,m+1)=LmTxz1*S(index((m+1):end)+m+1);

                                if m==1
                                    %Txz
                                    ATxz2{K+1}(:,m+1)=LmTxz2*deltaC(index((m+1):end)+m-1);
                                    BTxz2{K+1}(:,m+1)=LmTxz2*S(index((m+1):end)+m-1);

                                    %Tyz
                                    ATyz2{K+1}(:,m+1)=LmTyz2*deltaC(index((m+1):end)+m-1).*q;
                                    BTyz2{K+1}(:,m+1)=LmTyz2*S(index((m+1):end)+m-1).*q;
                                else
                                    %Txz
                                    ATxz2{K+1}(:,m+1)=0;
                                    BTxz2{K+1}(:,m+1)=0;

                                    %Tyz
                                    ATyz2{K+1}(:,m+1)=0;
                                    BTyz2{K+1}(:,m+1)=0;
                                end

                                %Txy
                                ATxy1{K+1}(:,m+1)=(LmTxy1*deltaC(index((m+1):end)+m+2)).*q;  
                                BTxy1{K+1}(:,m+1)=(LmTxy1*S(index((m+1):end)+m+2)).*q; 

                                ATxy3{K+1}(:,m+1)=0;
                                BTxy3{K+1}(:,m+1)=0;

                                %Tyz
                                ATyz1{K+1}(:,m+1)=(LmTyz1*deltaC(index((m+1):end)+m+1)).*q;  
                                BTyz1{K+1}(:,m+1)=(LmTyz1*S(index((m+1):end)+m+1)).*q; 
                            elseif m>nmax-2

                                if m==nmax
                                    %Txz
                                    ATxz1{K+1}(:,m+1)=0;
                                    BTxz1{K+1}(:,m+1)=0;

                                    %Tyz
                                    ATyz1{K+1}(:,m+1)=0;
                                    BTyz1{K+1}(:,m+1)=0;
                                else
                                    %Txz
                                    ATxz1{K+1}(:,m+1)=LmTxz1*deltaC(index((m+1):end)+m+1);
                                    BTxz1{K+1}(:,m+1)=LmTxz1*S(index((m+1):end)+m+1);

                                    %Tyz
                                    ATyz1{K+1}(:,m+1)=LmTyz1*deltaC(index((m+1):end)+m+1).*q;
                                    BTyz1{K+1}(:,m+1)=LmTyz1*S(index((m+1):end)+m+1).*q;
                                end

                                ATxz2{K+1}(:,m+1)=LmTxz2*deltaC(index((m+1):end)+m-1);
                                BTxz2{K+1}(:,m+1)=LmTxz2*S(index((m+1):end)+m-1);

                                %Txy
                                ATxy1{K+1}(:,m+1)=0;
                                BTxy1{K+1}(:,m+1)=0;

                                ATxy3{K+1}(:,m+1)=(LmTxy3*deltaC(index((m+1):end)+m-2)).*q;
                                BTxy3{K+1}(:,m+1)=(LmTxy3*S(index((m+1):end)+m-2)).*q;  

                                %Tyz
                                ATyz2{K+1}(:,m+1)=LmTyz2*deltaC(index((m+1):end)+m-1).*q;
                                BTyz2{K+1}(:,m+1)=LmTyz2*S(index((m+1):end)+m-1).*q;
                            else
                                ATxz1{K+1}(:,m+1)=LmTxz1*deltaC(index((m+1):end)+m+1);
                                BTxz1{K+1}(:,m+1)=LmTxz1*S(index((m+1):end)+m+1);

                                ATxz2{K+1}(:,m+1)=LmTxz2*deltaC(index((m+1):end)+m-1);
                                BTxz2{K+1}(:,m+1)=LmTxz2*S(index((m+1):end)+m-1);

                                %Txy
                                ATxy1{K+1}(:,m+1)=(LmTxy1*deltaC(index((m+1):end)+m+2)).*q;
                                BTxy1{K+1}(:,m+1)=(LmTxy1*S(index((m+1):end)+m+2)).*q;

                                ATxy3{K+1}(:,m+1)=(LmTxy3*deltaC(index((m+1):end)+m-2)).*q;
                                BTxy3{K+1}(:,m+1)=(LmTxy3*S(index((m+1):end)+m-2)).*q;

                                %Tyz
                                ATyz1{K+1}(:,m+1)=LmTyz1*deltaC(index((m+1):end)+m+1).*q;
                                BTyz1{K+1}(:,m+1)=LmTyz1*S(index((m+1):end)+m+1).*q;

                                ATyz2{K+1}(:,m+1)=LmTyz2*deltaC(index((m+1):end)+m-1).*q;
                                BTyz2{K+1}(:,m+1)=LmTyz2*S(index((m+1):end)+m-1).*q;
                            end
                        end
                        %Modified forward column method cannot be
                        %applied
                        
                    elseif volbapar(i)==10 %Geoid undulation

                        if m==nmax
                            amplH=zeros(length(fiG),nmax+1); %Damping factor
                            for n=0:nmax
                                amplH(:,n+1)=1./((R./r).^n);

                                % When computing H, there is no
                                % dumping factor (R./r).^n,
                                % therefore the matrix Pnm has to be
                                % devided by 1./((R./r).^n), since
                                % Pnm is the matrix of the MODIFIED
                                % fnALFS
                            end
                        end
                                
                        Lm=Pnm(:,1:(nmax-m+1));
                        LmH=Pnm(:,1:(nmax-m+1)).*amplH(:,(m+1):end);
                        
                        if volbaALFs==1 || volbaALFs==3
                            if m==nmax
                                AN1c=zeros(length(fiG),nmax+1);
                                BN1c=AN1c;
                                AH=AN1c;
                                BH=AN1c;
                            end  
                            
                            AN1c(:,m+1)=Lm*deltaCm;
                            BN1c(:,m+1)=Lm*Sm;  
                            AH(:,m+1)=LmH*HCm;
                            BH(:,m+1)=LmH*HSm;
                        elseif volbaALFs==2
                            if m==nmax
                                N1c=zeros(length(fiG),length(lambda));
                                H=N1c;
                            end

                            N1c=bsxfun(@times,N1c,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            H=bsxfun(@times,H,u)+(LmH*HCm*cos(m*lambda')+LmH*HSm*sin(m*lambda'));
                        end
                        
                    elseif volbapar(i)==11 %Gravitational potential

                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_V{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_V{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_V{K+1}=[ampl_V{K+1} sprintf('((0:nmax)+%d)',KK+1)];
                                    end
                                end
                                
                                ampl_V{K+1}=strrep(ampl_V{K+1},')(',').*(');
                                ampl_V{K+1}=eval(ampl_V{K+1});
                            end
                            
                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_V{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AV{K+1}=zeros(length(fiG),nmax+1);
                                    BV{K+1}=zeros(length(fiG),nmax+1);
                                end

                                AV{K+1}(:,m+1)=Lm*Cm;
                                BV{K+1}(:,m+1)=Lm*Sm;  
                            elseif volbaALFs==2
                                if m==nmax
                                    V{K+1}=zeros(length(fiG),length(lambda));
                                end

                                V{K+1}=bsxfun(@times,V{K+1},u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                        
                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Vrr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Vppll{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Vrr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Vppll{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Vrr{K+1}=[ampl_Vrr{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                        ampl_Vppll{K+1}=[ampl_Vppll{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Vrr{K+1}=strrep(ampl_Vrr{K+1},')(',').*(');
                                ampl_Vppll{K+1}=strrep(ampl_Vppll{K+1},')(',').*(');
                                ampl_Vrr{K+1}=eval(ampl_Vrr{K+1});
                                ampl_Vppll{K+1}=eval(ampl_Vppll{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vrr{K+1}(:,(m+1):end));
                            Lmff=bsxfun(@times,ddPnm(:,1:(nmax-m+1)),ampl_Vppll{K+1}(:,(m+1):end));
                            Lmll=bsxfun(@times,m^2*Pnm(:,1:(nmax-m+1)),ampl_Vppll{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AVrr{K+1}=zeros(length(fiG),nmax+1);
                                    BVrr{K+1}=zeros(length(fiG),nmax+1);
                                    AVff{K+1}=AVrr{K+1};
                                    BVff{K+1}=AVrr{K+1};
                                    AVll{K+1}=AVrr{K+1};
                                    BVll{K+1}=AVrr{K+1};
                                end

                                AVrr{K+1}(:,m+1)=Lm*Cm;
                                BVrr{K+1}(:,m+1)=Lm*Sm;  

                                AVff{K+1}(:,m+1)=Lmff*Cm;
                                BVff{K+1}(:,m+1)=Lmff*Sm;

                                AVll{K+1}(:,m+1)=Lmll*Cm;
                                BVll{K+1}(:,m+1)=Lmll*Sm;
                            elseif volbaALFs==2
                                if m==nmax
                                    Vrr{K+1}=zeros(length(fiG),length(lambda));
                                    Vff{K+1}=Vrr{K+1};
                                    Vll{K+1}=Vrr{K+1};
                                end

                                Vrr{K+1}=bsxfun(@times,Vrr{K+1},u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                Vff{K+1}=bsxfun(@times,Vff{K+1},u)+(Lmff*Cm*cos(m*lambda')+Lmff*Sm*sin(m*lambda'));
                                Vll{K+1}=bsxfun(@times,Vll{K+1},u)+(Lmll*Cm*cos(m*lambda')+Lmll*Sm*sin(m*lambda'));
                            end
                        end
                        
                    elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                        
                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Vrf{K+1}='((0:nmax)+1)';
                                    ampl_Vfl{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Vrf{K+1}='((0:nmax)+1)';
                                    ampl_Vfl{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Vrf{K+1}=[ampl_Vrf{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                        ampl_Vfl{K+1}=[ampl_Vfl{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Vrf{K+1}=strrep(ampl_Vrf{K+1},')(',').*(');
                                ampl_Vfl{K+1}=strrep(ampl_Vfl{K+1},')(',').*(');
                                ampl_Vrf{K+1}=eval(ampl_Vrf{K+1});
                                ampl_Vfl{K+1}=eval(ampl_Vfl{K+1});
                            end

                            Lmrf=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_Vrf{K+1}(:,(m+1):end));
                            Lmrl=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_Vrf{K+1}(:,(m+1):end));
                            Lmfl=bsxfun(@times,m*dPnm(:,1:(nmax-m+1)),ampl_Vfl{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AVrf{K+1}=zeros(length(fiG),nmax+1);
                                    BVrf{K+1}=zeros(length(fiG),nmax+1);
                                    AVrl{K+1}=AVrf{K+1};
                                    BVrl{K+1}=AVrf{K+1};
                                    AVfl{K+1}=AVrf{K+1};
                                    BVfl{K+1}=AVrf{K+1};
                                end

                                AVrf{K+1}(:,m+1)=Lmrf*Cm;
                                BVrf{K+1}(:,m+1)=Lmrf*Sm;  

                                AVrl{K+1}(:,m+1)=Lmrl*Cm;
                                BVrl{K+1}(:,m+1)=Lmrl*Sm;

                                AVfl{K+1}(:,m+1)=Lmfl*Cm;
                                BVfl{K+1}(:,m+1)=Lmfl*Sm;
                            elseif volbaALFs==2
                                if m==nmax
                                    Vrf{K+1}=zeros(length(fiG),length(lambda));
                                    Vrl{K+1}=Vrf{K+1};
                                    Vfl{K+1}=Vrf{K+1};
                                end

                                Vrf{K+1}=bsxfun(@times,Vrf{K+1},u)+(Lmrf*Cm*cos(m*lambda')+Lmrf*Sm*sin(m*lambda'));
                                Vrl{K+1}=bsxfun(@times,Vrl{K+1},u)+(Lmrl*Sm*cos(m*lambda')-Lmrl*Cm*sin(m*lambda'));
                                Vfl{K+1}=bsxfun(@times,Vfl{K+1},u)+(Lmfl*Sm*cos(m*lambda')-Lmfl*Cm*sin(m*lambda'));
                            end
                        end
                        
                    elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz

                        for K=0:TR
                            if m==nmax  
                                if K==0
                                    ampl_Vzz{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Vxxyy{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Vzz{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    ampl_Vxxyy{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Vzz{K+1}=[ampl_Vzz{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                        ampl_Vxxyy{K+1}=[ampl_Vxxyy{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end

                                ampl_Vzz{K+1}=strrep(ampl_Vzz{K+1},')(',').*(');
                                ampl_Vxxyy{K+1}=strrep(ampl_Vxxyy{K+1},')(',').*(');
                                ampl_Vzz{K+1}=eval(ampl_Vzz{K+1});
                                ampl_Vxxyy{K+1}=eval(ampl_Vxxyy{K+1});
                            end

                            if m<2
                                Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vzz{K+1}(:,(m+1):end));

                                %Vxx
                                bnm=((m:nmax)+m+1).*((m:nmax)+m+2)./2./(m+1);
                                Pnmxxyy=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vxxyy{K+1}(:,(m+1):end));
                                LmVxx1=bsxfun(@times,Pnmxxyy,bnm-((m:nmax)+1).*((m:nmax)+2));
                                LmVyy1=bsxfun(@times,Pnmxxyy,bnm);

                                if m==0
                                    anm=sqrt(2)/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                else
                                    anm=1./4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                end

                                LmVxx2=bsxfun(@times,Pnmxxyy,anm);
                                LmVxx3=zeros(length(fi),nmax-m+1); %Coefficients cnm are equal to zero, if m==0 a m==1
                            else
                                Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vzz{K+1}(:,(m+1):end));

                                %Vxx
                                bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                Pnmxxyy=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vxxyy{K+1}(:,(m+1):end));
                                LmVxx1=bsxfun(@times,Pnmxxyy,bnm-((m:nmax)+1).*((m:nmax)+2));
                                LmVyy1=bsxfun(@times,Pnmxxyy,bnm);

                                anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                LmVxx2=bsxfun(@times,Pnmxxyy,anm);

                                if m==2
                                    cnm=sqrt(2)/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                else
                                    cnm=1./4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                end

                                LmVxx3=bsxfun(@times,Pnmxxyy,cnm);                                    
                            end

                            if m==nmax
                                AVzz{K+1}=zeros(length(fiG),nmax+1);
                                BVzz{K+1}=zeros(length(fiG),nmax+1);
                                AVxx1{K+1}=AVzz{K+1};
                                BVxx1{K+1}=AVzz{K+1};
                                AVxx2{K+1}=AVzz{K+1};
                                BVxx2{K+1}=AVzz{K+1};
                                AVxx3{K+1}=AVzz{K+1};
                                BVxx3{K+1}=BVzz{K+1};
                                AVyy1{K+1}=AVzz{K+1};
                                BVyy1{K+1}=BVzz{K+1};
                            end

                            AVzz{K+1}(:,m+1)=Lm*Cm;
                            BVzz{K+1}(:,m+1)=Lm*Sm;  

                            %Vxx
                            AVxx1{K+1}(:,m+1)=LmVxx1*Cm;
                            BVxx1{K+1}(:,m+1)=LmVxx1*Sm;

                            if m<2
                                AVxx2{K+1}(:,m+1)=LmVxx2*C(index((m+1):end)+m+2);
                                BVxx2{K+1}(:,m+1)=LmVxx2*S(index((m+1):end)+m+2);

                                AVxx3{K+1}(:,m+1)=0;
                                BVxx3{K+1}(:,m+1)=0;
                            elseif m>nmax-2
                                AVxx2{K+1}(:,m+1)=0;
                                BVxx2{K+1}(:,m+1)=0;

                                AVxx3{K+1}(:,m+1)=LmVxx3*C(index((m+1):end)+m-2);
                                BVxx3{K+1}(:,m+1)=LmVxx3*S(index((m+1):end)+m-2);
                            else
                                AVxx2{K+1}(:,m+1)=LmVxx2*C(index((m+1):end)+m+2);
                                BVxx2{K+1}(:,m+1)=LmVxx2*S(index((m+1):end)+m+2);

                                AVxx3{K+1}(:,m+1)=LmVxx3*C(index((m+1):end)+m-2);
                                BVxx3{K+1}(:,m+1)=LmVxx3*S(index((m+1):end)+m-2);
                            end

                            %Vyy
                            AVyy1{K+1}(:,m+1)=LmVyy1*Cm;
                            BVyy1{K+1}(:,m+1)=LmVyy1*Sm;                                
                        end
                        %Modified forward column method cannot be
                        %applied
                        
                    elseif volbapar(i)==15 %Gravitational tensor Vxy_Vxz_Vyz
                        
                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Vnondiag{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_Vnondiag{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Vnondiag{K+1}=[ampl_Vnondiag{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Vnondiag{K+1}=strrep(ampl_Vnondiag{K+1},')(',').*(');
                                ampl_Vnondiag{K+1}=eval(ampl_Vnondiag{K+1});
                            end
                            
                            if m<2                                       
                                if m==0
                                    betanm=((m:nmax)+2)./2.*sqrt(1+ones(1,nmax+1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=zeros(1,nmax+1);

                                    %Vxy
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2); 
                                    
                                    Pnmnondiag=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ampl_Vnondiag{K+1}(:,(m+1):end));
                                    
                                    LmVxy1=bsxfun(@times,Pnmnondiag,dnm);
                                    LmVxy2=zeros(length(fi),nmax+1);
                                    LmVxy3=zeros(length(fi),nmax+1);

                                    %Vyz
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    LmVyz1=bsxfun(@times,Pnmnondiag,minm);

                                    LmVyz2=zeros(length(fi),nmax+1);
                                else
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2).*sqrt((m:nmax).*((m:nmax)+1)./2);

                                    Pnmnondiag=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ampl_Vnondiag{K+1}(:,(m+1):end));
                                    
                                    %Vxy
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);                                        
                                    LmVxy1=bsxfun(@times,Pnmnondiag,dnm);

                                    gnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+1).*sqrt((m:nmax)-1).*((m:nmax)+2);
                                    LmVxy2=bsxfun(@times,Pnmnondiag,gnm);

                                    LmVxy3=zeros(length(fi),nmax);

                                    %Vyz
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    LmVyz1=bsxfun(@times,Pnmnondiag,minm);

                                    LmVyz2=zeros(length(fi),nmax);
                                end

                                %Vxz
                                Pnmnondiag=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vnondiag{K+1}(:,(m+1):end));
                                LmVxz1=bsxfun(@times,Pnmnondiag,betanm);
                                LmVxz2=bsxfun(@times,Pnmnondiag,gamanm); 
                            else
                                %Vxz
                                betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);

                                Pnmnondiag=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vnondiag{K+1}(:,(m+1):end));
                                LmVxz1=bsxfun(@times,Pnmnondiag,betanm);
                                LmVxz2=bsxfun(@times,Pnmnondiag,gamanm);

                                %Vxy
                                dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                Pnmnondiag=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ampl_Vnondiag{K+1}(:,(m+1):end));
                                LmVxy1=bsxfun(@times,Pnmnondiag,dnm);

                                gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                LmVxy2=bsxfun(@times,Pnmnondiag,gnm);

                                if m==2
                                    %Vxy
                                    LmVxy3=zeros(length(fi),nmax-1);    
                                elseif m==3
                                    %Vxy
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                    LmVxy3=bsxfun(@times,Pnmnondiag,hnm);
                                else
                                    %Vxy
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                    LmVxy3=bsxfun(@times,Pnmnondiag,hnm);
                                end

                                %Vyz
                                minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                LmVyz1=bsxfun(@times,Pnmnondiag,minm);

                                ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                LmVyz2=bsxfun(@times,Pnmnondiag,ninm);
                            end

                            if m==nmax
                                AVxz1{K+1}=zeros(length(fiG),nmax+1);
                                BVxz1{K+1}=zeros(length(fiG),nmax+1);
                                AVxz2{K+1}=AVxz1{K+1};
                                BVxz2{K+1}=AVxz1{K+1};
                                AVxy1{K+1}=AVxz1{K+1};
                                BVxy1{K+1}=AVxz1{K+1};
                                AVxy2{K+1}=AVxz1{K+1};
                                BVxy2{K+1}=AVxz1{K+1};
                                AVxy3{K+1}=AVxz1{K+1};
                                BVxy3{K+1}=AVxz1{K+1};
                                AVyz1{K+1}=AVxz1{K+1};
                                BVyz1{K+1}=AVxz1{K+1};
                                AVyz2{K+1}=AVxz1{K+1};
                                BVyz2{K+1}=AVxz1{K+1};
                            end

                            AVxy2{K+1}(:,m+1)=(LmVxy2*Cm).*q;
                            BVxy2{K+1}(:,m+1)=(LmVxy2*Sm).*q;

                            if m<2
                                AVxz1{K+1}(:,m+1)=LmVxz1*C(index((m+1):end)+m+1);
                                BVxz1{K+1}(:,m+1)=LmVxz1*S(index((m+1):end)+m+1);

                                if m==1
                                    %Vxz
                                    AVxz2{K+1}(:,m+1)=LmVxz2*C(index((m+1):end)+m-1);
                                    BVxz2{K+1}(:,m+1)=LmVxz2*S(index((m+1):end)+m-1);

                                    %Vyz
                                    AVyz2{K+1}(:,m+1)=LmVyz2*C(index((m+1):end)+m-1).*q;
                                    BVyz2{K+1}(:,m+1)=LmVyz2*S(index((m+1):end)+m-1).*q;
                                else
                                    %Vxz
                                    AVxz2{K+1}(:,m+1)=0;
                                    BVxz2{K+1}(:,m+1)=0;

                                    %Vyz
                                    AVyz2{K+1}(:,m+1)=0;
                                    BVyz2{K+1}(:,m+1)=0;
                                end

                                %Vxy
                                AVxy1{K+1}(:,m+1)=(LmVxy1*C(index((m+1):end)+m+2)).*q;  
                                BVxy1{K+1}(:,m+1)=(LmVxy1*S(index((m+1):end)+m+2)).*q; 

                                AVxy3{K+1}(:,m+1)=0;
                                BVxy3{K+1}(:,m+1)=0;

                                %Vyz
                                AVyz1{K+1}(:,m+1)=(LmVyz1*C(index((m+1):end)+m+1)).*q;  
                                BVyz1{K+1}(:,m+1)=(LmVyz1*S(index((m+1):end)+m+1)).*q; 
                            elseif m>nmax-2

                                if m==nmax
                                    %Vxz
                                    AVxz1{K+1}(:,m+1)=0;
                                    BVxz1{K+1}(:,m+1)=0;

                                    %Vyz
                                    AVyz1{K+1}(:,m+1)=0;
                                    BVyz1{K+1}(:,m+1)=0;
                                else
                                    %Vxz
                                    AVxz1{K+1}(:,m+1)=LmVxz1*C(index((m+1):end)+m+1);
                                    BVxz1{K+1}(:,m+1)=LmVxz1*S(index((m+1):end)+m+1);

                                    %Vyz
                                    AVyz1{K+1}(:,m+1)=LmVyz1*C(index((m+1):end)+m+1).*q;
                                    BVyz1{K+1}(:,m+1)=LmVyz1*S(index((m+1):end)+m+1).*q;
                                end

                                AVxz2{K+1}(:,m+1)=LmVxz2*C(index((m+1):end)+m-1);
                                BVxz2{K+1}(:,m+1)=LmVxz2*S(index((m+1):end)+m-1);

                                %Vxy
                                AVxy1{K+1}(:,m+1)=0;
                                BVxy1{K+1}(:,m+1)=0;

                                AVxy3{K+1}(:,m+1)=(LmVxy3*C(index((m+1):end)+m-2)).*q;
                                BVxy3{K+1}(:,m+1)=(LmVxy3*S(index((m+1):end)+m-2)).*q;  

                                %Vyz
                                AVyz2{K+1}(:,m+1)=LmVyz2*C(index((m+1):end)+m-1).*q;
                                BVyz2{K+1}(:,m+1)=LmVyz2*S(index((m+1):end)+m-1).*q;
                            else
                                AVxz1{K+1}(:,m+1)=LmVxz1*C(index((m+1):end)+m+1);
                                BVxz1{K+1}(:,m+1)=LmVxz1*S(index((m+1):end)+m+1);

                                AVxz2{K+1}(:,m+1)=LmVxz2*C(index((m+1):end)+m-1);
                                BVxz2{K+1}(:,m+1)=LmVxz2*S(index((m+1):end)+m-1);

                                %Vxy
                                AVxy1{K+1}(:,m+1)=(LmVxy1*C(index((m+1):end)+m+2)).*q;
                                BVxy1{K+1}(:,m+1)=(LmVxy1*S(index((m+1):end)+m+2)).*q;

                                AVxy3{K+1}(:,m+1)=(LmVxy3*C(index((m+1):end)+m-2)).*q;
                                BVxy3{K+1}(:,m+1)=(LmVxy3*S(index((m+1):end)+m-2)).*q;

                                %Vyz
                                AVyz1{K+1}(:,m+1)=LmVyz1*C(index((m+1):end)+m+1).*q;
                                BVyz1{K+1}(:,m+1)=LmVyz1*S(index((m+1):end)+m+1).*q;

                                AVyz2{K+1}(:,m+1)=LmVyz2*C(index((m+1):end)+m-1).*q;
                                BVyz2{K+1}(:,m+1)=LmVyz2*S(index((m+1):end)+m-1).*q;
                            end
                        end
                        %Modified forward column method cannot be
                        %applied
                        
                    elseif volbapar(i)==16 %Gravity  

                        for K=0:TR
                            if m==nmax                                
                                if K==0
                                    ampl_Wr{K+1}='((0:nmax)+1)'; 
                                    ampl_g{K+1}='ones(1,nmax+1)'; 
                                else
                                    ampl_Wr{K+1}='((0:nmax)+1)';
                                    ampl_g{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Wr{K+1}=[ampl_Wr{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                        ampl_g{K+1}=[ampl_g{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end
                                end
                                
                                ampl_Wr{K+1}=strrep(ampl_Wr{K+1},')(',').*(');
                                ampl_g{K+1}=strrep(ampl_g{K+1},')(',').*(');
                                ampl_Wr{K+1}=eval(ampl_Wr{K+1});
                                ampl_g{K+1}=eval(ampl_g{K+1});
                            end

                            LmWr=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Wr{K+1}(:,(m+1):end));
                            LmWlambda=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_g{K+1}(:,(m+1):end));
                            LmWfi=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_g{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AWr{K+1}=zeros(length(fiG),nmax+1);
                                    BWr{K+1}=zeros(length(fiG),nmax+1);
                                    AWlambda{K+1}=AWr{K+1};
                                    BWlambda{K+1}=BWr{K+1};
                                    AWfi{K+1}=AWr{K+1};
                                    BWfi{K+1}=BWr{K+1};
                                end

                                AWr{K+1}(:,m+1)=LmWr*Cm;
                                BWr{K+1}(:,m+1)=LmWr*Sm;   
                                AWlambda{K+1}(:,m+1)=LmWlambda*Cm;
                                BWlambda{K+1}(:,m+1)=LmWlambda*Sm;
                                AWfi{K+1}(:,m+1)=LmWfi*Cm;
                                BWfi{K+1}(:,m+1)=LmWfi*Sm;
                            elseif volbaALFs==2
                                if m==nmax
                                    Wr{K+1}=zeros(length(fiG),length(lambda));
                                    Wlambda{K+1}=Wr{K+1};
                                    Wfi{K+1}=Wr{K+1};
                                end

                                Wr{K+1}=bsxfun(@times,Wr{K+1},u)+(LmWr*Cm*cos(m*lambda')+LmWr*Sm*sin(m*lambda'));
                                Wlambda{K+1}=bsxfun(@times,Wlambda{K+1},u)+(-LmWlambda*Cm*sin(m*lambda')+LmWlambda*Sm*cos(m*lambda'));
                                Wfi{K+1}=bsxfun(@times,Wfi{K+1},u)+(LmWfi*Cm*cos(m*lambda')+LmWfi*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==17 %Gravity sa

                        for K=0:TR
                            if m==nmax      
                                if K==0
                                    ampl_g_sa{K+1}='((0:nmax)+1)';
                                else
                                    ampl_g_sa{K+1}='((0:nmax)+1)';
                                    for KK=0:(K-1)
                                        ampl_g_sa{K+1}=[ampl_g_sa{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end
                                end
                                
                                ampl_g_sa{K+1}=strrep(ampl_g_sa{K+1},')(',').*(');
                                ampl_g_sa{K+1}=eval(ampl_g_sa{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_g_sa{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    Ag_sa{K+1}=zeros(length(fiG),nmax+1);
                                    Bg_sa{K+1}=zeros(length(fiG),nmax+1);
                                end

                                Ag_sa{K+1}(:,m+1)=Lm*Cm;
                                Bg_sa{K+1}(:,m+1)=Lm*Sm;                             
                            elseif volbaALFs==2
                                if m==nmax
                                   g_sa{K+1}=zeros(length(fiG),length(lambda)); 
                                end

                                g_sa{K+1}=bsxfun(@times,g_sa{K+1},u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==18 %Gravity potential                           

                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_W{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_W{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_W{K+1}=[ampl_W{K+1} sprintf('((0:nmax)+%d)',KK+1)];
                                    end
                                end
                                
                                ampl_W{K+1}=strrep(ampl_W{K+1},')(',').*(');
                                ampl_W{K+1}=eval(ampl_W{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_W{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AW{K+1}=zeros(length(fiG),nmax+1);
                                    BW{K+1}=zeros(length(fiG),nmax+1);
                                end

                                AW{K+1}(:,m+1)=Lm*Cm;
                                BW{K+1}(:,m+1)=Lm*Sm;  
                            elseif volbaALFs==2
                                if m==nmax
                                    W{K+1}=zeros(length(fiG),length(lambda));
                                end

                                W{K+1}=bsxfun(@times,W{K+1},u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==19 %Gravity anomaly sa

                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_anomalia_sa{K+1}='((0:nmax)-1)';    
                                else
                                    ampl_anomalia_sa{K+1}='((0:nmax)-1)';
                                    for KK=0:(K-1)
                                        ampl_anomalia_sa{K+1}=[ampl_anomalia_sa{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end
                                end
                                
                                ampl_anomalia_sa{K+1}=strrep(ampl_anomalia_sa{K+1},')(',').*(');
                                ampl_anomalia_sa{K+1}=eval(ampl_anomalia_sa{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_anomalia_sa{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    Aanomalia_sa{K+1}=zeros(length(fiG),nmax+1);
                                    Banomalia_sa{K+1}=zeros(length(fiG),nmax+1);
                                end

                                Aanomalia_sa{K+1}(:,m+1)=Lm*deltaCm;
                                Banomalia_sa{K+1}(:,m+1)=Lm*Sm;
                            elseif volbaALFs==2
                                if m==nmax
                                    anomalia_sa{K+1}=zeros(length(fiG),length(lambda));
                                end

                                anomalia_sa{K+1}=bsxfun(@times,anomalia_sa{K+1},u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==20 %Gravity disturbance

                        for K=0:TR
                            if m==nmax                                
                                if K==0
                                    ampl_Wrpor{K+1}='((0:nmax)+1)'; 
                                    ampl_gpor{K+1}='ones(1,nmax+1)'; 
                                else
                                    ampl_Wrpor{K+1}='((0:nmax)+1)';
                                    ampl_gpor{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_Wrpor{K+1}=[ampl_Wrpor{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                        ampl_gpor{K+1}=[ampl_gpor{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end
                                end
                                
                                ampl_Wrpor{K+1}=strrep(ampl_Wrpor{K+1},')(',').*(');
                                ampl_gpor{K+1}=strrep(ampl_gpor{K+1},')(',').*(');
                                ampl_Wrpor{K+1}=eval(ampl_Wrpor{K+1});
                                ampl_gpor{K+1}=eval(ampl_gpor{K+1});
                            end

                            LmWrpor=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Wrpor{K+1}(:,(m+1):end));
                            LmWlambdapor=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_gpor{K+1}(:,(m+1):end));
                            LmWfipor=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_gpor{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AWrpor{K+1}=zeros(length(fiG),nmax+1);
                                    BWrpor{K+1}=zeros(length(fiG),nmax+1);
                                    AWlambdapor{K+1}=AWrpor{K+1};
                                    BWlambdapor{K+1}=BWrpor{K+1};
                                    AWfipor{K+1}=AWrpor{K+1};
                                    BWfipor{K+1}=BWrpor{K+1};
                                end

                                AWrpor{K+1}(:,m+1)=LmWrpor*Cm;
                                BWrpor{K+1}(:,m+1)=LmWrpor*Sm;   
                                AWlambdapor{K+1}(:,m+1)=LmWlambdapor*Cm;
                                BWlambdapor{K+1}(:,m+1)=LmWlambdapor*Sm;
                                AWfipor{K+1}(:,m+1)=LmWfipor*Cm;
                                BWfipor{K+1}(:,m+1)=LmWfipor*Sm;

                                if m==0
                                    AUr{K+1}(:,m+1)=LmWrpor*CElm;
                                    AUfi{K+1}(:,m+1)=LmWfipor*CElm;
                                end
                            elseif volbaALFs==2
                                if m==nmax
                                    Wrpor{K+1}=zeros(length(fiG),length(lambda));
                                    Wlambdapor{K+1}=Wrpor{K+1};
                                    Wfipor{K+1}=Wrpor{K+1};
                                    Ur{K+1}=Wrpor{K+1};
                                    Ufi{K+1}=Wrpor{K+1};
                                end

                                Wrpor{K+1}=bsxfun(@times,Wrpor{K+1},u)+(LmWrpor*Cm*cos(m*lambda')+LmWrpor*Sm*sin(m*lambda'));
                                Wlambdapor{K+1}=bsxfun(@times,Wlambdapor{K+1},u)+(-LmWlambdapor*Cm*sin(m*lambda')+LmWlambdapor*Sm*cos(m*lambda'));
                                Wfipor{K+1}=bsxfun(@times,Wfipor{K+1},u)+(LmWfipor*Cm*cos(m*lambda')+LmWfipor*Sm*sin(m*lambda'));
                                
                                if m==0                                   
                                    Ur{K+1}=LmWrpor*CElm*cos(m*lambda');
                                    Ufi{K+1}=LmWfipor*CElm*cos(m*lambda');
                                end
                            end
                        end

                    elseif volbapar(i)==21 %Gravity disturbance sa

                        for K=0:TR
                            if m==nmax                                
                                if K==0
                                    ampl_porucha_sa{K+1}='((0:nmax)+1)';    
                                else
                                    ampl_porucha_sa{K+1}='((0:nmax)+1)';
                                    for KK=0:(K-1)
                                        ampl_porucha_sa{K+1}=[ampl_porucha_sa{K+1} sprintf('((0:nmax)+%d+1)',KK+1)];
                                    end
                                end
                                
                                ampl_porucha_sa{K+1}=strrep(ampl_porucha_sa{K+1},')(',').*(');
                                ampl_porucha_sa{K+1}=eval(ampl_porucha_sa{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_porucha_sa{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    Aporucha_sa{K+1}=zeros(length(fiG),nmax+1);
                                    Bporucha_sa{K+1}=zeros(length(fiG),nmax+1);
                                end

                                Aporucha_sa{K+1}(:,m+1)=Lm*deltaCm;
                                Bporucha_sa{K+1}(:,m+1)=Lm*Sm;
                            elseif volbaALFs==2
                                if m==nmax
                                    porucha_sa{K+1}=zeros(length(fiG),length(lambda));
                                end

                                porucha_sa{K+1}=bsxfun(@times,porucha_sa{K+1},u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==22 %Height anomaly ell

                        for K=0:TR
                            if m==nmax
                                if K==0
                                    ampl_zetaEll{K+1}='ones(1,nmax+1)';
                                else
                                    ampl_zetaEll{K+1}='ones(1,nmax+1)';
                                    for KK=0:(K-1)
                                        ampl_zetaEll{K+1}=[ampl_zetaEll{K+1} sprintf('((0:nmax)+%d)',KK+1)];
                                    end
                                end
                                
                                ampl_zetaEll{K+1}=strrep(ampl_zetaEll{K+1},')(',').*(');
                                ampl_zetaEll{K+1}=eval(ampl_zetaEll{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_zetaEll{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AzetaEl{K+1}=zeros(length(fiG),nmax+1);
                                    BzetaEl{K+1}=zeros(length(fiG),nmax+1);
                                end

                                AzetaEl{K+1}(:,m+1)=Lm*deltaCm;
                                BzetaEl{K+1}(:,m+1)=Lm*Sm;                            
                            elseif volbaALFs==2
                                if m==nmax
                                    zetaEl{K+1}=zeros(length(fiG),length(lambda));
                                end

                                zetaEl{K+1}=bsxfun(@times,zetaEl{K+1},u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==23 %Height anomaly
                        
                        if m==nmax
                            ampl_zeta_H=zeros(length(fiG),nmax+1); %Damping factor
                            for n=0:nmax
                                ampl_zeta_H(:,n+1)=1./((R./r).^n); 
                                % When computing H, there is no
                                % dumping factor (R./r).^n,
                                % therefore the matrix Pnm has to be
                                % devided by 1./((R./r).^n), since
                                % Pnm is the matrix of the MODIFIED
                                % fnALFS
                            end
                            
                            ampl_zeta_dg=((0:nmax)+1);
                        end

                        Lmdg=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_zeta_dg(:,(m+1):end));
                        Lm=Pnm(:,1:(nmax-m+1));
                        LmH=Pnm(:,1:(nmax-m+1)).*ampl_zeta_H(:,(m+1):end);
                        
                        if volbaALFs==1 || volbaALFs==3
                            if m==nmax
                                AN1c_zeta=zeros(length(fiG),nmax+1);
                                BN1c_zeta=AN1c_zeta;
                                Azetadg=AN1c_zeta;
                                Bzetadg=AN1c_zeta;
                                AH_zeta=AN1c_zeta;
                                BH_zeta=AN1c_zeta;
                                Azeta=AN1c_zeta;
                                Bzeta=AN1c_zeta;
                            end  
                            
                            AN1c_zeta(:,m+1)=Lm*deltaCm;
                            BN1c_zeta(:,m+1)=Lm*Sm;  
                            Azetadg(:,m+1)=Lmdg*deltaCm;
                            Bzetadg(:,m+1)=Lmdg*Sm;
                            AH_zeta(:,m+1)=LmH*HCm;
                            BH_zeta(:,m+1)=LmH*HSm;
                            Azeta(:,m+1)=Lm*deltaCm;
                            Bzeta(:,m+1)=Lm*Sm;
                        elseif volbaALFs==2
                            if m==nmax
                                zeta_N1c=zeros(length(fiG),length(lambda));
                                zeta_H=zeta_N1c;
                                zeta_dg=zeta_N1c;
                                zeta_zetaEl=zeta_N1c;
                            end

                            zeta_N1c=bsxfun(@times,zeta_N1c,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            zeta_H=bsxfun(@times,zeta_H,u)+(LmH*HCm*cos(m*lambda')+LmH*HSm*sin(m*lambda'));
                            zeta_dg=bsxfun(@times,zeta_dg,u)+(Lmdg*deltaCm*cos(m*lambda')+Lmdg*Sm*sin(m*lambda'));
                            zeta_zetaEl=bsxfun(@times,zeta_zetaEl,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                        end
                                                        
                    elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                        for K=0:TR
                            if m==nmax       
                                if K==0
                                    ampl_T_rr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                else
                                    ampl_T_rr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    for KK=0:(K-1)
                                        ampl_T_rr{K+1}=[ampl_T_rr{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_T_rr{K+1}=strrep(ampl_T_rr{K+1},')(',').*(');
                                ampl_T_rr{K+1}=eval(ampl_T_rr{K+1});
                            end
                            
                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_T_rr{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AT_rr{K+1}=zeros(length(fiG),nmax+1);
                                    BT_rr{K+1}=zeros(length(fiG),nmax+1);
                                end

                                AT_rr{K+1}(:,m+1)=Lm*deltaCm;
                                BT_rr{K+1}(:,m+1)=Lm*Sm;                             
                            elseif volbaALFs==2
                                if m==nmax
                                    T_rr{K+1}=zeros(length(fiG),length(lambda));
                                end

                                T_rr{K+1}=bsxfun(@times,T_rr{K+1},u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end

                    elseif volbapar(i)==25 %Second radial derivative of disturbing potential

                        for K=0:TR
                            if m==nmax                               
                                if K==0
                                    ampl_Wrr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                else
                                    ampl_Wrr{K+1}='((0:nmax)+1).*((0:nmax)+2)';
                                    for KK=0:(K-1)
                                        ampl_Wrr{K+1}=[ampl_Wrr{K+1} sprintf('((0:nmax)+%d+2)',KK+1)];
                                    end
                                end
                                
                                ampl_Wrr{K+1}=strrep(ampl_Wrr{K+1},')(',').*(');
                                ampl_Wrr{K+1}=eval(ampl_Wrr{K+1});
                            end

                            Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Wrr{K+1}(:,(m+1):end));

                            if volbaALFs==1 || volbaALFs==3
                                if m==nmax
                                    AWrr{K+1}=zeros(length(fiG),nmax+1);
                                    BWrr{K+1}=zeros(length(fiG),nmax+1);
                                end

                                AWrr{K+1}(:,m+1)=Lm*Cm;
                                BWrr{K+1}(:,m+1)=Lm*Sm;                           
                            elseif volbaALFs==2
                                if m==nmax
                                    Wrr{K+1}=zeros(length(fiG),length(lambda));
                                end

                                Wrr{K+1}=bsxfun(@times,Wrr{K+1},u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                            end
                        end
                    end
                end                                       
            end                

            clear Lm dLm Pnm dPnm ddPnm Cm Sm C CEl CElm deltaC ...
                deltaCm S u t q q2 index tu qu enm singdPnm singddPnm
            
            clear ampl_eta ampl_ksi ampl_Teta ampl_Tksi ampl_T ampl_Trr ...
                ampl_Tppll ampl_Trf ampl_Tfl ampl_Tzz ampl_Txxyy ampl_Tnondiag ...
                amplH ampl_V ampl_Vrr ampl_Vppll ampl_Vrf ampl_Vfl ampl_Vzz ...
                ampl_Vxxyy ampl_Vnondiag ampl_Wr ampl_g ampl_g_sa ...
                ampl_W ampl_anomalia_sa ampl_Wrpor ampl_gpor ampl_porucha_sa ...
                ampl_zetaEll ampl_H ampl_dg ampl_T_rr ampl_Wrr
            
            if volbaALFs==1 || volbaALFs==3
                cosla=cos((0:nmax)'*lambda');
                sinla=sin((0:nmax)'*lambda');
            end
            
            %Update of the progress bar
            set(progressbar,'string',...
                'Progress: Loading DTM file...',...
                'fontsize',8); drawnow;
         
            if any(volbapar==10) || any(volbapar==23)
            else
                %Loading DMR (txt or mat file)                
                if strcmp(loadname(end-3:end),'.mat')
                    DTM=load([loadadresar,loadname]);
                    DTM=struct2cell(DTM);
                    DTM=cell2mat(DTM);
                else
                    DTM=load([loadadresar,loadname]);
                end

                if load_matrix==1
                    [rows_DTM,cols_DTM]=size(DTM);
                    if length(fi)*length(lambda)==rows_DTM*cols_DTM
                        DTM=DTM(end:-1:1,:);
                    else
                        errordlg(sprintf('The number of points in the input file containing heights of the irregular surface (%d) is not equal to the number of points in the entered grid (%d).',rows_DTM*cols_DTM,length(fi)*length(lambda)),...
                                'Error in irregular surface selection');
                        error('The number of points in the input file containing heights of the irregular surface is not equal to the number of points in the entered grid.')
                    end
                elseif load_vector==1
                    [rows_DTM,cols_DTM]=size(DTM);
                    
                    if cols_DTM~=1
                        errordlg('The input file containing heights of the irregular surface does not have the required structure of a column vector.',...
                            'Error in irregular surface selection')
                        error('The input file containing heights of the irregular surface does not have the required structure of a column vector.')
                    end
                    
                    %Reshape column vector to matrix
                    if length(fi)*length(lambda)==length(DTM)
                        DTM=reshape(DTM,length(fi),length(lambda));
                    else
                        errordlg(sprintf('The number of points in the input file containing heights of the irregular surface (%d) is not equal to the number of points in the entered grid (%d).',rows_DTM*cols_DTM,length(fi)*length(lambda)),...
                                'Error in irregular surface selection');
                        error('The number of points in the input file containing heights of the irregular surface is not equal to the number of points in the entered grid.')
                    end
                end 
                
                %Mean value of input heights
                meanDTM=mean(mean(DTM));
            end 
            
            if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==10) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20) || any(volbapar==22) || any(volbapar==23) || any(volbapar==25)
                %Computation of spherical radii and geocentric
                %latitudes of the points at irregular surface. They are
                %stored and restored when computing e.g. centrifugal
                %potential, etc.
                if coord==1
                    [LAMBDA,FIG]=meshgrid(lambda,fiG); 
                    [X,Y,Z]=sph2cart(LAMBDA,FIG,DTM); %DTM contains spherical radii
                elseif coord==0
                    [LAMBDA,FI]=meshgrid(lambda,fi); 
                    
                    if any(volbapar==10) || any(volbapar==23)
                        [X,Y,Z]=geodetic2ecef(FI,LAMBDA,0,[aEl eEl]); %Ellipsoidal height si in the case of geoid/height anomaly set to zero
                    else
                        [X,Y,Z]=geodetic2ecef(FI,LAMBDA,DTM,[aEl eEl]); %DTM contains ellipsoidal heights
                    end
                end                                                               
                clear FI FIG LAMBDA      

                if any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20)
                    rDTM=sqrt(X.^2+Y.^2+Z.^2);
                end
                
                if coord==0
                    if any(volbapar==2) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20) || any(volbapar==25)
                        fiGDTM=atan(Z./sqrt(X.^2+Y.^2));
                    end
                end

                %Computation of the normal gravity for eta, xi, Theta,
                %Geoid undulation, zeta el, zeta (at the irregular surface
                %except for the geoid undulation and height anomaly)
                if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==10) || any(volbapar==22) || any(volbapar==23)
                    %Update of the progress bar
                    set(progressbar,'string',...
                        'Progress: Comp. normal gravity...',...
                        'fontsize',8); drawnow;
            
                    bEl=aEl*sqrt(1-eEl^2);
                    EEl=sqrt(aEl^2-bEl^2);                   

                    %Computation of ellipsoidal harmonic coordinates
                    ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                    betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                    clear X Y Z

                    wgama=sqrt((ugama.^2+EEl^2*sin(betagama).^2)./(ugama.^2+EEl^2));
                    qgama=1/2*((1+(3*ugama.^2)./EEl^2).*atan(EEl./ugama)-3*ugama./EEl);
                    qgama_=3*(1+(ugama.^2)./EEl^2).*(1-ugama./EEl.*atan(EEl./ugama))-1;
                    qgama0=1/2*((1+(3*bEl^2)/EEl^2)*atan(EEl/bEl)-3*bEl/EEl);

                    gamau=-1./wgama.*(GMEl./(ugama.^2+EEl.^2)+(omegaEl.^2.*aEl.^2.*EEl)./(ugama.^2+EEl.^2).*(qgama_./qgama0).*(1./2.*sin(betagama).^2-1./6)-omegaEl.^2.*ugama.*cos(betagama).^2);
                    gamabeta=-1./wgama.*(-(omegaEl.^2.*aEl.^2.*qgama)./(sqrt(ugama.^2+EEl.^2).*qgama0)+omegaEl.^2.*sqrt(ugama.^2+EEl.^2)).*sin(betagama).*cos(betagama);

                    clear ugama betagama wgama qgama qgama_ qgama0

                    gamaP=sqrt(gamau.^2+gamabeta.^2);

                    clear gamau gamabeta
                end
                clear X Y Z
            end
            
            %Update of the progress bar
            set(progressbar,'string',...
                'Progress: Matrix multiplications...',...
                'fontsize',8); drawnow;
            
            %% Final computation of functionals of the geopotential             
            for i=1:pocetpar
                if volbapar(i)==1                
                elseif volbapar(i)==2 %Deflection of the vertical eta
                    
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            eta=-Aeta{K+1}*sinla+Beta{K+1}*cosla;                                
                            
                            Aeta{K+1}=[];
                            Beta{K+1}=[];
                            
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),eta).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),eta).*bsxfun(@minus,DTM,r).^K;
                            end                                                              
                        elseif volbaALFs==2
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),eta{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),eta{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            eta{K+1}=[];
                        end
                    end

                    clear Aeta Beta eta
                    Pg(fi==pi/2 | fi==-pi/2,:)=0; %#ok<*AGROW>
                    if coord==1
                        Pg=Pg./bsxfun(@times,gamaP,cos(fiG))*(180/pi)*3600;
                        Pg=Pg(:);
                    elseif coord==0
                        Pg=Pg(:)./(gamaP(:).*cos(fiGDTM(:)))*(180/pi)*3600;
                    end
                    
                elseif volbapar(i)==3 %Deflection of the vertical xi
                    
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            ksi=Aksi{K+1}*cosla+Bksi{K+1}*sinla;
                            
                            Aksi{K+1}=[]; 
                            Bksi{K+1}=[];
                            
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),ksi).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),ksi).*bsxfun(@minus,DTM,r).^K;
                            end
                        elseif volbaALFs==2
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),ksi{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),ksi{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            ksi{K+1}=[];
                        end
                    end
                    
                    clear Aksi Bksi ksi                        
                    Pg=Pg(:)./(gamaP(:))*(180/pi)*3600;
                    
                elseif volbapar(i)==4 %Deflection of the vertical Theta
                    
                    Thetaeta=0;
                    Thetaksi=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            Teta=-ATeta{K+1}*sinla+BTeta{K+1}*cosla;
                            Tksi=ATksi{K+1}*cosla+BTksi{K+1}*sinla;
                            
                            ATksi{K+1}=[]; 
                            BTksi{K+1}=[];
                            ATeta{K+1}=[]; 
                            BTeta{K+1}=[];
                            
                            if coord==0
                                Thetaksi=Thetaksi+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Tksi).*bsxfun(@minus,DTM,h).^K;
                                Thetaeta=Thetaeta+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Teta).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Thetaksi=Thetaksi+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Tksi).*bsxfun(@minus,DTM,r).^K;
                                Thetaeta=Thetaeta+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Teta).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                        elseif volbaALFs==2
                            if coord==0
                                Thetaeta=Thetaeta+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Teta{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                                Thetaksi=Thetaksi+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Tksi{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Thetaeta=Thetaeta+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Teta{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                                Thetaksi=Thetaksi+1/factorial(K)*bsxfun(@times,(-1)^(K+1)*GM./r.^(K+2),Tksi{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            Teta{K+1}=[];
                            Tksi{K+1}=[];
                        end
                    end
                    
                    clear ATeta BTeta ATksi BTksi Teta Tksi
                    
                    Thetaeta(fi==pi/2 | fi==-pi/2,:)=0;
                    if coord==1
                        Thetaeta=Thetaeta./bsxfun(@times,gamaP,cos(fiG))*(180/pi)*3600;
                        Thetaeta=Thetaeta(:);
                    elseif coord==0
                        Thetaeta=Thetaeta(:)./(gamaP(:).*cos(fiGDTM(:)))*(180/pi)*3600;                
                    end
                    Thetaksi=Thetaksi(:)./(gamaP(:))*(180/pi)*3600;

                    Thetaalfa=atan2(Thetaeta,Thetaksi);
                    Thetaalfa(Thetaalfa<0)=Thetaalfa(Thetaalfa<0)+2*pi;
   
                    Pg=[sqrt(Thetaeta.^2+Thetaksi.^2) 180/pi*(Thetaalfa)];

                    clear Thetaeta Thetaksi Thetaalfa
                    
                elseif volbapar(i)==5 %Disturbing potential

                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            T=AT{K+1}*cosla+BT{K+1}*sinla;
                            
                            AT{K+1}=[];
                            BT{K+1}=[];
                            
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),T).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),T).*bsxfun(@minus,DTM,r).^K;
                            end                                
                        elseif volbaALFs==2
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),T{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                if nmin0==1 %Zero degree term
                                    Pg=Pg+bsxfun(@plus,(-1)^K*(GM-GMEl)./r.^(K+1),1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),T{K+1}*1e280)).*bsxfun(@minus,DTM,r).^K;
                                else
                                    Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),T{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                                end
                            end
                            
                            T{K+1}=[];
                        end
                    end
                   
                    clear AT BT T
                    Pg=Pg(:);   
                    
                elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                    
                    clear Lmff Lmll                            
                    
                    Trr_tenzor=0;
                    Tff_tenzor=0;
                    Tll_tenzor=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            %Trr                                                                
                            if coord==0
                                Trr_tenzor=Trr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATrr{K+1}*cosla+BTrr{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Trr_tenzor=Trr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATrr{K+1}*cosla+BTrr{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end                                
                            
                            ATrr{K+1}=[];
                            BTrr{K+1}=[];
                             
                            %Tff
                            if coord==0
                                Tff_tenzor=Tff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATff{K+1}*cosla+BTff{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tff_tenzor=Tff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATff{K+1}*cosla+BTff{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
                                                            
                            ATff{K+1}=[];
                            BTff{K+1}=[];
                            
                            %Tll
                            ATll{K+1}(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                            BTll{K+1}(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                            
                            if coord==0
                                Tll_tenzor=Tll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATll{K+1}*cosla+BTll{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tll_tenzor=Tll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATll{K+1}*cosla+BTll{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            ATll{K+1}=[];
                            BTll{K+1}=[];
                        elseif volbaALFs==2
                            %Trr                                                                
                            if coord==0
                                Trr_tenzor=Trr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Trr{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Trr_tenzor=Trr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Trr{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            Trr{K+1}=[];
                            
                            %Tff
                            if coord==0
                                Tff_tenzor=Tff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Tff{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tff_tenzor=Tff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Tff{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Tff{K+1}=[];
                            
                            %Tll
                            Tll{K+1}(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                            
                            if coord==0
                                Tll_tenzor=Tll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Tll{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tll_tenzor=Tll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Tll{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Tll{K+1}=[];
                        end
                    end
                    
                    clear ATrr BTrr ATff BTff ATll BTll
                    
                    Pg=Trr_tenzor(:)*10^9;
                    clear Trr_tenzor
                    Pg=[Pg Tff_tenzor(:)*10^9];
                    clear Tff_tenzor
                    if coord==1
                        Tll_tenzor=bsxfun(@rdivide,Tll_tenzor,cos(fiG).^2);
                        Pg=[Pg -Tll_tenzor(:)*10^9];
                    elseif coord==0
                        Pg=[Pg (-Tll_tenzor(:)./cos(fiGDTM(:)).^2)*10^9];
                    end
                    clear Tll_tenzor 

                elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                    
                    clear Lmrf Lmrl Lmfl
                    
                    Trf_tenzor=0;
                    Trl_tenzor=0;
                    Tfl_tenzor=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            %Trf                               
                            if coord==0
                                Trf_tenzor=Trf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATrf{K+1}*cosla+BTrf{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Trf_tenzor=Trf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATrf{K+1}*cosla+BTrf{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            ATrf{K+1}=[];
                            BTrf{K+1}=[];
                            
                            %Trl
                            if coord==0
                                Trl_tenzor=Trl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATrl{K+1}*sinla+BTrl{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Trl_tenzor=Trl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATrl{K+1}*sinla+BTrl{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            ATrl{K+1}=[];
                            BTrl{K+1}=[];
                            
                            %Tfl
                            ATfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            BTfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            
                            if coord==0
                                Tfl_tenzor=Tfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATfl{K+1}*sinla+BTfl{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tfl_tenzor=Tfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATfl{K+1}*sinla+BTfl{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            ATfl{K+1}=[];
                            BTfl{K+1}=[];
                        elseif volbaALFs==2                               
                            %Trf                               
                            if coord==0
                                Trf_tenzor=Trf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Trf{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Trf_tenzor=Trf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Trf{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            Trf{K+1}=[];
                            
                            %Trl
                            if coord==0
                                Trl_tenzor=Trl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Trl{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Trl_tenzor=Trl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Trl{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Trl{K+1}=[];

                            %Tfl
                            ATfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            BTfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            
                            if coord==0
                                Tfl_tenzor=Tfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Tfl{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tfl_tenzor=Tfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Tfl{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Tfl{K+1}=[];
                        end
                    end
                    
                    clear ATrf BTrf ATrl BTrl ATfl BTfl
                    
                    Pg=-Trf_tenzor(:)*10^9;
                    clear Trf_tenzor
                    if coord==1
                        Trl_tenzor=bsxfun(@rdivide,Trl_tenzor,cos(fiG));
                        Pg=[Pg -Trl_tenzor(:)*10^9];
                        clear Trl_tenzor
                        Tfl_tenzor=bsxfun(@rdivide,Tfl_tenzor,cos(fiG));
                        Pg=[Pg Tfl_tenzor(:)*10^9];
                        clear Tfl_tenzor
                    elseif coord==0
                        Pg=[Pg -Trl_tenzor(:)./cos(fiGDTM(:))*10^9];
                        clear Trl_tenzor
                        Pg=[Pg Tfl_tenzor(:)./cos(fiGDTM(:))*10^9];
                        clear Tfl_tenzor
                    end
                    
                elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                    
                    clear anm bnm cnm LmTxx1 LmTxx2 LmTxx3 ...
                        LmTyy1 ampl_Tzz Pnmxxyy
                    
                    Tzz_tenzor=0;
                    Txx1_tenzor=0;
                    Txx2_tenzor=0;
                    Tyy1_tenzor=0;
                    Txx3_tenzor=0;
                    
                    GMindT=zeros(K+1,1);
                    GMindT(1)=1;
                    for K=1:TR
                        GMindT(K+1)=GMindT(K)+(K+1);
                    end
                    
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3                                     
                            %Tzz                                                               
                            if coord==0
                                Tzz_tenzor=Tzz_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATzz{K+1}*cosla+BTzz{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tzz_tenzor=Tzz_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATzz{K+1}*cosla+BTzz{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end                                
                            
                            ATzz{K+1}=[];
                            BTzz{K+1}=[];

                            %Txx
                            if coord==0
                                Txx1_tenzor=Txx1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxx1{K+1}*cosla+BTxx1{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txx1_tenzor=Txx1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxx1{K+1}*cosla+BTxx1{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            ATxx1{K+1}=[];
                            BTxx1{K+1}=[];

                            if coord==0
                                Txx2_tenzor=Txx2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxx2{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BTxx2{K+1}*[sinla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txx2_tenzor=Txx2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxx2{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BTxx2{K+1}*[sinla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            ATxx2{K+1}=[];
                            BTxx2{K+1}=[];

                            Txx_tenzor=Txx1_tenzor+Txx2_tenzor;

                            %Tyy
                            if coord==0
                                Tyy1_tenzor=Tyy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATyy1{K+1}*cosla+BTyy1{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tyy1_tenzor=Tyy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATyy1{K+1}*cosla+BTyy1{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            ATyy1{K+1}=[];
                            BTyy1{K+1}=[];

                            Tyy_tenzor=Tyy1_tenzor+Txx2_tenzor;

                            if coord==0
                                Txx3_tenzor=Txx3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxx3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BTxx3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txx3_tenzor=Txx3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxx3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BTxx3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            ATxx3{K+1}=[];
                            BTxx3{K+1}=[];

                            Txx_tenzor=Txx_tenzor+Txx3_tenzor;                                                                                                                       
                            Tyy_tenzor=Tyy_tenzor+Txx3_tenzor;

                        elseif volbaALFs==2
                        end
                    end
                    
                    clear GMindT ATzz BTzz ATxx1 BTxx1 ATxx2 BTxx2 ATyy1 ...
                        BTyy1 ATxx3 BTxx3 Txx1_tenzor Txx2_tenzor ...
                        Tyy1_tenzor Txx3_tenzor
                    
                    Pg=Txx_tenzor(:)*10^9;
                    clear Txx_tenzor
                    Pg=[Pg -Tyy_tenzor(:)*10^9];
                    clear Tyy_tenzor
                    Pg=[Pg Tzz_tenzor(:)*10^9];
                    clear Tzz_tenzor

                elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                    
                    clear gamanm betanm gnm hnm dnm minm ninm ...
                        LmTxz1 LmTxz2 LmTxy1 LmTxy2 LmTxy3 ...
                        LmTyz1 LmTyz2 Pnmnondiag
                    
                    Txz1_tenzor=0;
                    Txz2_tenzor=0;
                    Txy1_tenzor=0;
                    Txy2_tenzor=0;
                    Txy3_tenzor=0;
                    Tyz1_tenzor=0;
                    Tyz2_tenzor=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3

                            %Txz                               
                            if coord==0
                                Txz1_tenzor=Txz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BTxz1{K+1}*[sinla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txz1_tenzor=Txz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BTxz1{K+1}*[sinla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            ATxz1{K+1}=[];
                            BTxz1{K+1}=[];

                            if coord==0
                                Txz2_tenzor=Txz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BTxz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txz2_tenzor=Txz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),ATxz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BTxz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            ATxz2{K+1}=[];
                            BTxz2{K+1}=[];

                            Txz_tenzor=Txz1_tenzor+Txz2_tenzor;

                            %Txy
                            if coord==0
                                Txy1_tenzor=Txy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATxy1{K+1}*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BTxy1{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txy1_tenzor=Txy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATxy1{K+1}*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BTxy1{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            ATxy1{K+1}=[];
                            BTxy1{K+1}=[];

                            if coord==0
                                Txy2_tenzor=Txy2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATxy2{K+1}*sinla+BTxy2{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txy2_tenzor=Txy2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATxy2{K+1}*sinla+BTxy2{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end
  
                            ATxy2{K+1}=[];
                            BTxy2{K+1}=[];

                            if coord==0
                                Txy3_tenzor=Txy3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATxy3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BTxy3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Txy3_tenzor=Txy3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATxy3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BTxy3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            ATxy3{K+1}=[];
                            BTxy3{K+1}=[];

                            Txy_tenzor=Txy1_tenzor+Txy2_tenzor+Txy3_tenzor;

                            %Tyz
                            if coord==0
                                Tyz1_tenzor=Tyz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATyz1{K+1}*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BTyz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tyz1_tenzor=Tyz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATyz1{K+1}*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BTyz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end

                            ATyz1{K+1}=[];
                            BTyz1{K+1}=[];

                            if coord==0
                                Tyz2_tenzor=Tyz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATyz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BTyz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Tyz2_tenzor=Tyz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-ATyz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BTyz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]).*bsxfun(@minus,DTM,r).^K;
                            end

                            ATyz2{K+1}=[];
                            BTyz2{K+1}=[];

                            Tyz_tenzor=Tyz1_tenzor+Tyz2_tenzor;
                        elseif volbaALFs==2
                        end
                    end
                    
                    clear ATxz1 BTxz1 ATxz2 BTxz2 ATxy1 BTxy1 ATxy2 ...
                        BTxy2 ATxy3 BTxy3 ATyz1 BTyz1 ATyz2 BTyz2 ...
                        Txz1_tenzor Txz2_tenzor Txy1_tenzor Txy2_tenzor ...
                        Txy3_tenzor Tyz1_tenzor Tyz2_tenzor                       
                     
                    Pg=Txy_tenzor(:)*10^9;
                    clear Txy_tenzor
                    Pg=[Pg Txz_tenzor(:)*10^9];
                    clear Txz_tenzor
                    Pg=[Pg Tyz_tenzor(:)*10^9];
                    clear Tyz_tenzor

                elseif volbapar(i)==10 %Geoid undulation
                    
                    clear HC HCm HS HSm LmH
                    
                    if volbaALFs==1 || volbaALFs==3
                        H=AH*cosla+BH*sinla;
                        clear AH BH
                        N1c=AN1c*cosla+BN1c*sinla;
                        clear AN1c BN1c
                    elseif volbaALFs==2
                        H=H*1e280;
                        N1c=N1c*1e280;
                    end                            
                    
                    N1c=bsxfun(@times,GM./r,N1c)./gamaP;
                    H(H<0)=H(H<0)*0; %H is set to zero in the areas of oceans and seas
                   
                    G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                    ro=2670; %Density of the crust

                    Pg=N1c-bsxfun(@times,(2*pi*G*ro*H.^2),1./gamaP);
                    Pg=Pg(:);
                    
                    clear H N1c
                elseif volbapar(i)==11 %Gravitational potential
                    
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            V=AV{K+1}*cosla+BV{K+1}*sinla;
                            
                            AV{K+1}=[];
                            BV{K+1}=[];
                            
                            if coord==0                                    
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),V).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),V).*bsxfun(@minus,DTM,r).^K;
                            end                                                                                              
                        elseif volbaALFs==2
                            if coord==0                                    
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),V{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),V{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            V{K+1}=[];
                        end                            
                    end
                    
                    clear V AV BV
                    Pg=Pg(:);
                    
                elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                    
                    clear Lmff Lmll                            

                    Vrr_tenzor=0;
                    Vff_tenzor=0;
                    Vll_tenzor=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            %Vrr                                                                
                            if coord==0
                                Vrr_tenzor=Vrr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVrr{K+1}*cosla+BVrr{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vrr_tenzor=Vrr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVrr{K+1}*cosla+BVrr{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end                                
                            
                            AVrr{K+1}=[];
                            BVrr{K+1}=[];
                             
                            %Vff
                            if coord==0
                                Vff_tenzor=Vff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVff{K+1}*cosla+BVff{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vff_tenzor=Vff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVff{K+1}*cosla+BVff{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
                                                            
                            AVff{K+1}=[];
                            BVff{K+1}=[];
                            
                            %Vll
                            AVll{K+1}(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                            BVll{K+1}(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                            
                            if coord==0
                                Vll_tenzor=Vll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVll{K+1}*cosla+BVll{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vll_tenzor=Vll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVll{K+1}*cosla+BVll{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            AVll{K+1}=[];
                            BVll{K+1}=[];
                        elseif volbaALFs==2
                            %Vrr                                                                
                            if coord==0
                                Vrr_tenzor=Vrr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vrr{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vrr_tenzor=Vrr_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vrr{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            Vrr{K+1}=[];
                            
                            %Vff
                            if coord==0
                                Vff_tenzor=Vff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vff{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vff_tenzor=Vff_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vff{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Vff{K+1}=[];
                            
                            %Vll
                            Vll{K+1}(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                            
                            if coord==0
                                Vll_tenzor=Vll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vll{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vll_tenzor=Vll_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vll{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Vll{K+1}=[];
                        end
                    end
                    
                    clear AVrr BVrr AVff BVff AVll BVll Vrr Vff Vll
                    
                    Pg=Vrr_tenzor(:)*10^9;
                    clear Vrr_tenzor
                    Pg=[Pg Vff_tenzor(:)*10^9];
                    clear Vff_tenzor
                    if coord==1
                        Vll_tenzor=bsxfun(@rdivide,Vll_tenzor,cos(fiG).^2);
                        Pg=[Pg -Vll_tenzor(:)*10^9];
                    elseif coord==0
                        Pg=[Pg (-Vll_tenzor(:)./cos(fiGDTM(:)).^2)*10^9];
                    end
                    clear Vll_tenzor  
    
                elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                    
                    clear Lmrf Lmrl Lmfl

                    Vrf_tenzor=0;
                    Vrl_tenzor=0;
                    Vfl_tenzor=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            %Vrf                               
                            if coord==0
                                Vrf_tenzor=Vrf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVrf{K+1}*cosla+BVrf{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vrf_tenzor=Vrf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVrf{K+1}*cosla+BVrf{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            AVrf{K+1}=[];
                            BVrf{K+1}=[];
                            
                            %Vrl
                            if coord==0
                                Vrl_tenzor=Vrl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVrl{K+1}*sinla+BVrl{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vrl_tenzor=Vrl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVrl{K+1}*sinla+BVrl{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end  
                            
                            AVrl{K+1}=[];
                            BVrl{K+1}=[];
                            
                            %Vfl
                            AVfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            BVfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            
                            if coord==0
                                Vfl_tenzor=Vfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVfl{K+1}*sinla+BVfl{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vfl_tenzor=Vfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVfl{K+1}*sinla+BVfl{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            AVfl{K+1}=[];
                            BVfl{K+1}=[];
                        elseif volbaALFs==2                               
                            %Vrf                               
                            if coord==0
                                Vrf_tenzor=Vrf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vrf{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vrf_tenzor=Vrf_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vrf{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            Vrf{K+1}=[];
                            
                            %Vrl
                            if coord==0
                                Vrl_tenzor=Vrl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vrl{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vrl_tenzor=Vrl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vrl{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Vrl{K+1}=[];

                            %Vfl
                            AVfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            BVfl{K+1}(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                            
                            if coord==0
                                Vfl_tenzor=Vfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vfl{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vfl_tenzor=Vfl_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Vfl{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Vfl{K+1}=[];
                        end
                    end
                    
                    clear AVrf BVrf AVrl BVrl AVfl BVfl Vrf Vrl Vfl
                    
                    Pg=-Vrf_tenzor(:)*10^9;
                    clear Vrf_tenzor
                    if coord==1
                        Vrl_tenzor=bsxfun(@rdivide,Vrl_tenzor,cos(fiG));
                        Pg=[Pg -Vrl_tenzor(:)*10^9];
                        clear Vrl_tenzor
                        Vfl_tenzor=bsxfun(@rdivide,Vfl_tenzor,cos(fiG));
                        Pg=[Pg Vfl_tenzor(:)*10^9];
                        clear Vfl_tenzor
                    elseif coord==0
                        Pg=[Pg -Vrl_tenzor(:)./cos(fiGDTM(:))*10^9];
                        clear Vrl_tenzor
                        Pg=[Pg Vfl_tenzor(:)./cos(fiGDTM(:))*10^9];
                        clear Vfl_tenzor
                    end
                    
                elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                    
                    clear anm bnm cnm LmVxx1 LmVxx2 LmVxx3 ...
                        LmVyy1 ampl_Vzz Pnmxxyy
                    
                    Vzz_tenzor=0;
                    Vxx1_tenzor=0;
                    Vxx2_tenzor=0;
                    Vyy1_tenzor=0;
                    Vxx3_tenzor=0;
                    
                    GMind=zeros(K+1,1);
                    GMind(1)=1;
                    for K=1:TR
                        GMind(K+1)=GMind(K)+(K+1);
                    end
            
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            %Vzz                                                               
                            if coord==0
                                Vzz_tenzor=Vzz_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVzz{K+1}*cosla+BVzz{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vzz_tenzor=Vzz_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVzz{K+1}*cosla+BVzz{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end                                
                            
                            AVzz{K+1}=[];
                            BVzz{K+1}=[];

                            %Vxx
                            if coord==0
                                Vxx1_tenzor=Vxx1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxx1{K+1}*cosla+BVxx1{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxx1_tenzor=Vxx1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxx1{K+1}*cosla+BVxx1{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end 
        
                            AVxx1{K+1}=[];
                            BVxx1{K+1}=[];

                            if coord==0
                                Vxx2_tenzor=Vxx2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxx2{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BVxx2{K+1}*[sinla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxx2_tenzor=Vxx2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxx2{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BVxx2{K+1}*[sinla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            AVxx2{K+1}=[];
                            BVxx2{K+1}=[];

                            Vxx_tenzor=Vxx1_tenzor+Vxx2_tenzor;

                            %Tyy
                            if coord==0
                                Vyy1_tenzor=Vyy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVyy1{K+1}*cosla+BVyy1{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vyy1_tenzor=Vyy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVyy1{K+1}*cosla+BVyy1{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end  
     
                            AVyy1{K+1}=[];
                            BVyy1{K+1}=[];

                            Vyy_tenzor=Vyy1_tenzor+Vxx2_tenzor;

                            if coord==0
                                Vxx3_tenzor=Vxx3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxx3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BVxx3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxx3_tenzor=Vxx3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxx3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BVxx3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            AVxx3{K+1}=[];
                            BVxx3{K+1}=[];

                            Vxx_tenzor=Vxx_tenzor+Vxx3_tenzor;                                                                                                                       
                            Vyy_tenzor=Vyy_tenzor+Vxx3_tenzor;

                        elseif volbaALFs==2
                        end
                    end
                    
                    clear GMind AVzz BVzz AVxx1 BVxx1 AVxx2 BVxx2 ...
                        AVyy1 BVyy1 AVxx3 BVxx3 Vxx1_tenzor ...
                        Vxx2_tenzor Vyy1_tenzor Vxx3_tenzor
                    
                    Pg=Vxx_tenzor(:)*10^9;
                    clear Vxx_tenzor
                    Pg=[Pg -Vyy_tenzor(:)*10^9];
                    clear Vyy_tenzor
                    Pg=[Pg Vzz_tenzor(:)*10^9];
                    clear Vzz_tenzor
                    
                elseif volbapar(i)==15 %Gravitational tensor Vxy_Vxz_Vyz
                    
                    clear gamanm betanm gnm hnm dnm minm ninm ...
                        LmVxz1 LmVxz2 LmVxy1 LmVxy2 LmVxy3 ...
                        LmVyz1 LmVyz2 Pnmnondiag
                    
                    Vxz1_tenzor=0;
                    Vxz2_tenzor=0;
                    Vxy1_tenzor=0;
                    Vxy2_tenzor=0;
                    Vxy3_tenzor=0;
                    Vyz1_tenzor=0;
                    Vyz2_tenzor=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3

                            %Vxz                               
                            if coord==0
                                Vxz1_tenzor=Vxz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BVxz1{K+1}*[sinla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxz1_tenzor=Vxz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BVxz1{K+1}*[sinla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            AVxz1{K+1}=[];
                            BVxz1{K+1}=[];

                            if coord==0
                                Vxz2_tenzor=Vxz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BVxz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxz2_tenzor=Vxz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),AVxz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BVxz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            AVxz2{K+1}=[];
                            BVxz2{K+1}=[];

                            Vxz_tenzor=Vxz1_tenzor+Vxz2_tenzor;

                            %Vxy
                            if coord==0
                                Vxy1_tenzor=Vxy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVxy1{K+1}*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BVxy1{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxy1_tenzor=Vxy1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVxy1{K+1}*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BVxy1{K+1}*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            AVxy1{K+1}=[];
                            BVxy1{K+1}=[];

                            if coord==0
                                Vxy2_tenzor=Vxy2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVxy2{K+1}*sinla+BVxy2{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxy2_tenzor=Vxy2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVxy2{K+1}*sinla+BVxy2{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end
  
                            AVxy2{K+1}=[];
                            BVxy2{K+1}=[];

                            if coord==0
                                Vxy3_tenzor=Vxy3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVxy3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BVxy3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vxy3_tenzor=Vxy3_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVxy3{K+1}*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BVxy3{K+1}*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            AVxy3{K+1}=[];
                            BVxy3{K+1}=[];

                            Vxy_tenzor=Vxy1_tenzor+Vxy2_tenzor+Vxy3_tenzor;

                            %Vyz
                            if coord==0
                                Vyz1_tenzor=Vyz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVyz1{K+1}*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BVyz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vyz1_tenzor=Vyz1_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVyz1{K+1}*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BVyz1{K+1}*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]).*bsxfun(@minus,DTM,r).^K;
                            end

                            AVyz1{K+1}=[];
                            BVyz1{K+1}=[];

                            if coord==0
                                Vyz2_tenzor=Vyz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVyz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BVyz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Vyz2_tenzor=Vyz2_tenzor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),-AVyz2{K+1}*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BVyz2{K+1}*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]).*bsxfun(@minus,DTM,r).^K;
                            end

                            AVyz2{K+1}=[];
                            BVyz2{K+1}=[];

                            Vyz_tenzor=Vyz1_tenzor+Vyz2_tenzor;
                        elseif volbaALFs==2
                        end
                    end
                    
                    clear AVxz1 BVxz1 AVxz2 BVxz2 AVxy1 BVxy1 AVxy2 ...
                        BVxy2 AVxy3 BVxy3 AVyz1 BVyz1 AVyz2 BVyz2 ...
                        Vxz1_tenzor Vxz2_tenzor Vxy1_tenzor Vxy2_tenzor ...
                        Vxy3_tenzor Vyz1_tenzor Vyz2_tenzor                       
                     
                    Pg=Vxy_tenzor(:)*10^9;
                    clear Vxy_tenzor
                    Pg=[Pg Vxz_tenzor(:)*10^9];
                    clear Vxz_tenzor
                    Pg=[Pg Vyz_tenzor(:)*10^9];
                    clear Vyz_tenzor
                    
                elseif volbapar(i)==16 %Gravity
                    clear LmWr LmWlambda LmWfi
                    
                    gWr=0;
                    gWlambda=0;
                    gWfi=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            %Wr
                            if coord==0                                    
                                gWr=gWr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWr{K+1}*cosla+BWr{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWr=gWr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWr{K+1}*cosla+BWr{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            AWr{K+1}=[];
                            BWr{K+1}=[];
                            
                            %Wlambda
                            if coord==0
                                gWlambda=gWlambda+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),-AWlambda{K+1}*sinla+BWlambda{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWlambda=gWlambda+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),-AWlambda{K+1}*sinla+BWlambda{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            AWlambda{K+1}=[];
                            BWlambda{K+1}=[];
                            
                            %Wfi
                            if coord==0
                                gWfi=gWfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWfi{K+1}*cosla+BWfi{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWfi=gWfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWfi{K+1}*cosla+BWfi{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            AWfi{K+1}=[];
                            BWfi{K+1}=[];
                        elseif volbaALFs==2
                            %Wr
                            if coord==0                                    
                                gWr=gWr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wr{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWr=gWr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wr{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 

                            Wr{K+1}=[];
                            
                            %Wlambda
                            if coord==0
                                gWlambda=gWlambda+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wlambda{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWlambda=gWlambda+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wlambda{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Wlambda{K+1}=[];
                            
                            %Wfi
                            if coord==0
                                gWfi=gWfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wfi{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWfi=gWfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wfi{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 
                            
                            Wfi{K+1}=[];
                        end                           
                    end

                    clear AWr BWr AWfi BWfi AWlambda BWlambda Wr Wfi Wlambda
                    
                    if coord==1
                        gWr=-gWr+omegaEl^2.*bsxfun(@times,rDTM,cos(fiG).^2);
                        gWr=gWr(:);
                        gWlambda=bsxfun(@rdivide,gWlambda,cos(fiG));
                        gWlambda=gWlambda(:);
                        gWfi=gWfi-omegaEl^2.*bsxfun(@times,rDTM,cos(fiG).*sin(fiG));
                        gWfi=gWfi(:);
                    elseif coord==0
                        gWr=-gWr(:)+omegaEl^2.*rDTM(:).*cos(fiGDTM(:)).^2;
                        gWlambda=gWlambda(:)./cos(fiGDTM(:));
                        gWfi=gWfi(:)-omegaEl^2.*rDTM(:).*cos(fiGDTM(:)).*sin(fiGDTM(:));
                    end
                    
                    Pg=sqrt(gWr.^2+gWlambda.^2+gWfi.^2)*10^5;

                    clear gWr gWlambda gWfi
                elseif volbapar(i)==17 %Gravity sa
                    
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            g_sa=Ag_sa{K+1}*cosla+Bg_sa{K+1}*sinla;
                            
                            Ag_sa{K+1}=[];
                            Bg_sa{K+1}=[];
                            
                            if coord==0                                    
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),g_sa).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),g_sa).*bsxfun(@minus,DTM,r).^K;
                            end                               
                        elseif volbaALFs==2
                            if coord==0                                    
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),g_sa{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),g_sa{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            g_sa{K+1}=[];
                        end                 
                    end
                    
                    clear Ag_sa Bg_sa g_sa
                    if coord==1
                        Pg=sqrt((-Pg+omegaEl^2.*bsxfun(@times,rDTM,cos(fiG).^2)).^2)*10^5;
                        Pg=Pg(:);
                    elseif coord==0
                        Pg=sqrt((-Pg(:)+omegaEl^2.*rDTM(:).*cos(fiGDTM(:)).^2).^2)*10^5;
                    end
                    
                elseif volbapar(i)==18 %Gravity potential
                    
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            W=AW{K+1}*cosla+BW{K+1}*sinla;
                            
                            AW{K+1}=[];
                            BW{K+1}=[];
                            
                            if coord==0  
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),W).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),W).*bsxfun(@minus,DTM,r).^K;
                            end                                                                
                        elseif volbaALFs==2
                            if coord==0                                    
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),W{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),W{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            W{K+1}=[];
                        end                            
                    end
                    
                    clear W AW BW
                    if coord==1
                        Pg=Pg+1/2*omegaEl.^2.*bsxfun(@times,rDTM.^2,cos(fiG).^2);
                        Pg=Pg(:);
                    elseif coord==0
                        Pg=Pg(:)+1/2*omegaEl.^2.*rDTM(:).^2.*cos(fiGDTM(:)).^2;
                    end

                elseif volbapar(i)==19 %Gravity anomaly sa
                    
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            anomalia_sa=Aanomalia_sa{K+1}*cosla+Banomalia_sa{K+1}*sinla;

                            Aanomalia_sa{K+1}=[];
                            Banomalia_sa{K+1}=[];
                            
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),anomalia_sa).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),anomalia_sa).*bsxfun(@minus,DTM,r).^K;
                            end                                                                
                        elseif volbaALFs==2                               
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),anomalia_sa{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),anomalia_sa{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            anomalia_sa{K+1}=[];
                        end

                    end
                    
                    clear Aanomalia_sa Banomalia_sa anomalia_sa 
                    Pg=Pg(:)*10^5;
                    
                elseif volbapar(i)==20 %Gravity disturbance
                    clear LmWrpor LmWlambdapor LmWfipor
                    
                    gWrpor=0;
                    gWfipor=0;
                    gWlambdapor=0;
                    gUr=0;
                    gUfi=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            %Wr
                            if coord==0                                    
                                gWrpor=gWrpor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWrpor{K+1}*cosla+BWrpor{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWrpor=gWrpor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWrpor{K+1}*cosla+BWrpor{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end 

                            AWrpor{K+1}=[];
                            BWrpor{K+1}=[];

                            %Wlambda
                            if coord==0
                               gWlambdapor=gWlambdapor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),-AWlambdapor{K+1}*sinla+BWlambdapor{K+1}*cosla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWlambdapor=gWlambdapor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),-AWlambdapor{K+1}*sinla+BWlambdapor{K+1}*cosla).*bsxfun(@minus,DTM,r).^K;
                            end 

                            AWlambdapor{K+1}=[];
                            BWlambdapor{K+1}=[];

                            %Wfi
                            if coord==0
                                gWfipor=gWfipor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWfipor{K+1}*cosla+BWfipor{K+1}*sinla).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWfipor=gWfipor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AWfipor{K+1}*cosla+BWfipor{K+1}*sinla).*bsxfun(@minus,DTM,r).^K;
                            end 

                            AWfipor{K+1}=[];
                            BWfipor{K+1}=[];
                            
                            %Ur
                            if coord==0    
                                gUr=gUr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AUr{K+1}*cos(0*lambda')).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1                                
                                gUr=gUr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AUr{K+1}*cos(0*lambda')).*bsxfun(@minus,DTM,r).^K;
                            end 

                            AUr{K+1}=[];
                            
                            %Ufi
                            if coord==0
                                gUfi=gUfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AUfi{K+1}*cos(0*lambda')).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gUfi=gUfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),AUfi{K+1}*cos(0*lambda')).*bsxfun(@minus,DTM,r).^K;
                            end 

                            AUfi{K+1}=[];
                        elseif volbaALFs==2
                            %Wr
                            if coord==0
                                gWrpor=gWrpor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wrpor{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1                                
                                gWrpor=gWrpor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wrpor{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 


                            Wrpor{K+1}=[];
                            
                            %Wlambda
                            if coord==0
                               gWlambdapor=gWlambdapor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wlambdapor{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWlambdapor=gWlambdapor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wlambdapor{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 

                            Wlambdapor{K+1}=[];

                            %Wfi
                            if coord==0
                                gWfipor=gWfipor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wfipor{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gWfipor=gWfipor+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Wfipor{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 

                            Wfipor{K+1}=[];
                            
                            %Ur
                            if coord==0
                                gUr=gUr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Ur{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1                                
                                gUr=gUr+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Ur{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 

                            Ur{K+1}=[];
                            
                            %Ufi
                            if coord==0
                                gUfi=gUfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Ufi{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                gUfi=gUfi+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),Ufi{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end 

                            Ufi{K+1}=[];
                        end
                    end

                    clear AWrpor BWrpor AWfipor BWfipor AWlambdapor ...
                        BWlambdapor Wrpor Wfipor Wlambdapor Wrpor ...
                        Wlambdapor Wfipor AUr AUfi Ufi Ur           
                    
                    %gWrpor
                    if coord==1
                        gWrpor=sqrt(bsxfun(@plus,-gWrpor,omegaEl^2.*r.*(cos(fiG).^2)).^2);
                    elseif coord==0
                        gWrpor=sqrt((-gWrpor+omegaEl^2.*rDTM.*(cos(fiGDTM).^2)).^2);
                    end
                    gWrpor=gWrpor(:);

                    %gWlambdapor
                    if coord==1
                        gWlambdapor=bsxfun(@rdivide,gWlambdapor,cos(fiG));
                        gWlambdapor=gWlambdapor(:);
                    elseif coord==0
                        gWlambdapor=gWlambdapor./cos(fiGDTM);                     
                    end
                    gWlambdapor=gWlambdapor(:);
             
                    %gWfipor
                    if coord==1
                        gWfipor=bsxfun(@plus,gWfipor,-omegaEl^2.*r.*cos(fiG).*sin(fiG));
                    elseif coord==0
                        gWfipor=gWfipor-omegaEl^2.*rDTM.*cos(fiGDTM).*sin(fiGDTM);
                    end                    
                    gWfipor=gWfipor(:);

                    %gUr
                    if coord==1
                        gUr=sqrt(bsxfun(@plus,-gUr,omegaEl^2.*r.*(cos(fiG).^2)).^2);
                    elseif coord==0
                        gUr=sqrt((-gUr+omegaEl^2.*rDTM.*(cos(fiGDTM).^2)).^2);
                    end
                    gUr=gUr(:);
                    
                    %gUfi
                    if coord==1
                        gUfi=bsxfun(@plus,gUfi,-omegaEl^2.*r.*cos(fiG).*sin(fiG));
                    elseif coord==0
                        gUfi=gUfi-omegaEl^2.*rDTM.*cos(fiGDTM).*sin(fiGDTM);
                    end                    
                    gUfi=gUfi(:);
                   
                    Pg=(sqrt(gWrpor.^2+gWlambdapor.^2+gWfipor.^2)-sqrt(gUr.^2+gUfi.^2))*10^5;
    
                    clear gWrpor gWlambdapor gWfipor gUfi gUr                         
                elseif volbapar(i)==21 %Gravity disturbance sa
                
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            porucha_sa=Aporucha_sa{K+1}*cosla+Bporucha_sa{K+1}*sinla; 
                            
                            Aporucha_sa{K+1}=[];
                            Bporucha_sa{K+1}=[];
                             
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),porucha_sa).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),porucha_sa).*bsxfun(@minus,DTM,r).^K;
                            end                                                              
                        elseif volbaALFs==2                                
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),porucha_sa{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+2),porucha_sa{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                             
                            porucha_sa{K+1}=[];
                        end                     
                    end
                     
                    clear Aporucha_sa Bporucha_sa porucha_sa 
                    Pg=Pg(:)*10^5;
                     
                elseif volbapar(i)==22 %Height anomaly ell

                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            zetaEl=AzetaEl{K+1}*cosla+BzetaEl{K+1}*sinla;
                             
                            AzetaEl{K+1}=[];
                            BzetaEl{K+1}=[];
                             
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),zetaEl).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),zetaEl).*bsxfun(@minus,DTM,r).^K;
                            end                                                             
                        elseif volbaALFs==2
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),zetaEl{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+1),zetaEl{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                             
                            zetaEl{K+1}=[];
                        end
                    end
                     
                    clear zetaEl AzetaEl BzetaEl
                    Pg=Pg(:)./gamaP(:);

                elseif volbapar(i)==23 %Height anomaly
                     
                    clear HC HCm HS HSm LmH Lmdg
                     
                    if volbaALFs==1 || volbaALFs==3
                        zeta_H=AH_zeta*cosla+BH_zeta*sinla;
                        clear AH_zeta BH_zeta
                        zeta_N1c=AN1c_zeta*cosla+BN1c_zeta*sinla;
                        clear AN1c_zeta BN1c_zeta
                        zeta_dg=Azetadg*cosla+Bzetadg*sinla;
                        clear Azetadg Bzetadg
                        zeta_zetaEl=Azeta*cosla+Bzeta*sinla;
                        clear Azeta Bzeta  
                    elseif volbaALFs==2
                        zeta_H=zeta_H*1e280;
                        zeta_N1c=zeta_N1c*1e280;
                        zeta_dg=zeta_dg*1e280;
                        zeta_zetaEl=zeta_zetaEl*1e280;
                    end                                 
                     
                    zeta_N1c=bsxfun(@times,GM./r,zeta_N1c)./gamaP;
                    zeta_H(zeta_H<0)=zeta_H(zeta_H<0)*0; %H is set to zero in the areas of oceans and seas   
                    zeta_zetaEl=bsxfun(@times,GM./r,zeta_zetaEl)./gamaP;                    
                    zeta_dg=bsxfun(@times,GM./r.^2,zeta_dg);
                     
                    G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                    ro=2670; %Density of the crust
                     
                    zeta_N=(zeta_N1c-bsxfun(@times,(2*pi*G*ro*zeta_H.^2),1./gamaP));

                    Pg=zeta_zetaEl-bsxfun(@times,zeta_dg.*(zeta_H+zeta_N),1./gamaP);
                    Pg=Pg(:);
                     
                    clear zeta_N1c zeta_H zeta_zetaEl zeta_dg zeta_N
                elseif volbapar(i)==24 %Second radial derivative of disturbing potential
                     
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            T_rr=AT_rr{K+1}*cosla+BT_rr{K+1}*sinla;
                            
                            AT_rr{K+1}=[];
                            BT_rr{K+1}=[];
                            
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),T_rr).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),T_rr).*bsxfun(@minus,DTM,r).^K;
                            end
                        elseif volbaALFs==2
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),T_rr{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),T_rr{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                             
                            T_rr{K+1}=[];
                        end
                    end
                    
                    clear T_rr AT_rr BT_rr
                    Pg=Pg(:)*10^9;

                elseif volbapar(i)==25 %Second radial derivative of gravity potential
                     
                    Pg=0;
                    for K=TR:-1:0
                        if volbaALFs==1 || volbaALFs==3
                            Wrr=AWrr{K+1}*cosla+BWrr{K+1}*sinla;
                            
                            AWrr{K+1}=[];
                            BWrr{K+1}=[];
                            
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Wrr).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Wrr).*bsxfun(@minus,DTM,r).^K;
                            end                                
                        elseif volbaALFs==2
                            if coord==0
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Wrr{K+1}*1e280).*bsxfun(@minus,DTM,h).^K;
                            elseif coord==1
                                Pg=Pg+1/factorial(K)*bsxfun(@times,(-1)^K*GM./r.^(K+3),Wrr{K+1}*1e280).*bsxfun(@minus,DTM,r).^K;
                            end
                            
                            Wrr{K+1}=[];
                        end
                    end
                     
                    clear AWrr BWrr Wrr
                    if coord==1
                        Pg=bsxfun(@plus,Pg,omegaEl^2.*cos(fiG).^2);
                        Pg=Pg(:)*10^9;
                    elseif coord==0
                        Pg=(Pg(:)+omegaEl^2.*cos(fiGDTM(:)).^2)*10^9;
                    end
                end
    
                if i==1
                    P=Pg;
                    clear Pg
                else
                    P=[P Pg];
                    if i==pocetpar
                        clear Pg
                    end
                end
            end

            clear sinla cosla r gamaP rDTM fiGDTM 
            
            %Points with no heights should be indicated by -9999 or NaN
            if geoid~=1
                P(DTM(:)==-9999,:)=-9999;
            end
    
            %Update of the progress bar
            set(progressbar,'string','','fontsize',8); drawnow;
                        
            set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;               
                
            cas=toc; %Stop clock

            %% Export data
            %==============================================================
                        
            if coord==1 %Spherical coordinates
                fi=fiG; %Ellipsoidal latitude is replaced by the spherical one
                clear fiG %Redundand vector fiG is deleted
            elseif coord==0 %Ellipsoidal coordinates
                clear fiG %Spherical latitude is deleted, since exported is
                %to be the ellipsoidal one
            end
            
            if get(findobj('tag','export'),'value')==1 || get(findobj('tag','datamat'),'value')==1
                   
                set(findobj('tag','hlasky'),'string',...
                    'Creating data file...','fontsize',8,...
                    'foregroundcolor','k'); drawnow;
                
                [lambda,fi]=meshgrid(180/pi*(lambda),180/pi*(fi));
                [j_fi,i_fi]=size(fi);
                fi=fi(:);
                lambda=lambda(:);                                        

                try
                    if geoid==0
                        output=[fi';lambda';DTM(:)';P'];
                    elseif geoid==1
                        output=[fi';lambda';P'];
                    end
                    clear DTM
                    
                    if get(findobj('tag','export'),'value')==1 %Export to data file
                        [rows_output, cols_output]=size(output);

                        exp1=fopen([outadresar,[outname '.txt']],'w');
                        fprintf(exp1,...
                            [repmat('% 0.12e ',1,rows_output-1) '% 0.12e\n'],output);                   
                        fclose(exp1);
                        set(findobj('tag','hlasky'),'string',...
                        '','fontsize',8,...
                        'foregroundcolor','k'); drawnow;  
                    end
    
                    if get(findobj('tag','datamat'),'value')==1 %Export to .mat file
                        output=output'; %#ok<*NASGU>
                        save([outadresar,[outname '.mat']],'output');
                    end
                
                    clear output
                    
                    export_error=0;
                 catch err
                     export_error=1;
                     
                     save([outadresar,[outname '_phi.mat']],'fi');
                     save([outadresar,[outname '_lambda.mat']],'lambda');
                     save([outadresar,[outname '_functionals.mat']],'P');
                 end
            end
            
            %% Export report
            %============================================================== 
            
            if get(findobj('tag','report'),'value')==1
               
                set(findobj('tag','hlasky'),'string',...
                    'Creating report file...','fontsize',8,...
                    'foregroundcolor','k'); drawnow;
                
                date_time=clock;
                date_time=datestr(date_time);
                
                reportname=[outname '_Report' '.txt'];
                exp=fopen([outadresar,reportname],'w');
                                
                if get(findobj('tag','export'),'value')==0 && get(findobj('tag','datamat'),'value')==0
                    fi=180/pi*(fi);
                    lambda=180/pi*(lambda);
                    j_fi=length(fi);
                    i_fi=length(lambda);                    
                end  

                fprintf(exp,'Software                                           \tisGrafLab 1.1.2\n');
                fprintf(exp,'Generating date                                    \t');
                fprintf(exp,'%s',date_time(1:11));
                fprintf(exp,'\nGenerating time                                  \t');
                fprintf(exp,'%s',date_time(13:20));
                
                fprintf(exp,'\nGeopotential model file                          \t');
                fprintf(exp,'%s',[GGMadresar,GGMname]);
                
                if get(findobj('tag','P1'),'value')==10 || get(findobj('tag','P2'),'value')==10 || get(findobj('tag','P3'),'value')==10 || get(findobj('tag','P4'),'value')==10 || get(findobj('tag','P1'),'value')==23 || get(findobj('tag','P2'),'value')==23 || get(findobj('tag','P3'),'value')==23 || get(findobj('tag','P4'),'value')==23
                    fprintf(exp,'\nDigital terrain model file                      \t');
                    fprintf(exp,'%s',[loadadresarDMR,loadnameDMR]);
                    fprintf(exp,'\nNewtonian gravitational constant (m^3.kg^-1.s^-2)  \t');
                    fprintf(exp,'%1.5e',G);
                    fprintf(exp,'\nDensity of the crust (kg.m^-3)                     \t');
                    fprintf(exp,'%4.0f',ro);
                else
                    fprintf(exp,'\nInput file of irregular surface heights            \t');
                    fprintf(exp,'%s',[loadadresar,loadname]);
                end

                fprintf(exp,'\nGM of the geopotential model (m^3.s^-2)          \t');                
                fprintf(exp,'%1.9e',GM);
                fprintf(exp,'\nR of the geopotential model (m)                  \t');
                fprintf(exp,'%1.9e',R);
                fprintf(exp,'\nMinimum used degree                              \t');
                fprintf(exp,'%0.0f',nmin);
                fprintf(exp,'\nMaximum used degree                              \t');
                fprintf(exp,'%0.0f',nmax);  
                fprintf(exp,'\nReference ellipsoid                              \t');

                if get(findobj('tag','ell'),'value') == 1;
                    fprintf(exp,'GRS80');
                elseif get(findobj('tag','ell'),'value') == 2;
                    fprintf(exp,'WGS84');
                end
                                
                fprintf(exp,'\nType of the input coordinates                    \t');
                if coord==1
                    fprintf(exp,'Spherical');
                elseif coord==0
                    fprintf(exp,'Ellipsoidal');
                end
                fprintf(exp,'\nLatitude limit North (deg)                       \t');
                fprintf(exp,'%0.9f',max(fi));
                fprintf(exp,'\nLatitude limit South (deg)                       \t');
                fprintf(exp,'%0.9f',min(fi));
                fprintf(exp,'\nLongitude limit West (deg)                       \t');
                fprintf(exp,'%0.9f',min(lambda));
                fprintf(exp,'\nLongitude limit East (deg)                       \t');
                fprintf(exp,'%0.9f',max(lambda));

                fprintf(exp,'\nLatitude parallels                               \t');
                fprintf(exp,'%0.0f',j_fi);
                fprintf(exp,'\nLongitude parallels                              \t');
                fprintf(exp,'%0.0f',i_fi);
                fprintf(exp,'\nNumber of grid points                            \t');
                fprintf(exp,'%0.0f',j_fi*i_fi);                                     
                fprintf(exp,'\nGrid height above the reference surface (m)      \t');
                    
                if coord==1 %Spherical coordinates
                    fprintf(exp,'%0.3f',hsph);
                elseif coord==0 %Ellipsoidal coordinates
                    fprintf(exp,'%0.3f',h(1));
                end
                clear h
                
                if get(findobj('tag','P1'),'value')==10 || get(findobj('tag','P2'),'value')==10 || get(findobj('tag','P3'),'value')==10 || get(findobj('tag','P4'),'value')==10 || get(findobj('tag','P1'),'value')==23 || get(findobj('tag','P2'),'value')==23 || get(findobj('tag','P3'),'value')==23 || get(findobj('tag','P4'),'value')==23
                else
                    fprintf(exp,'\nMean height (m)                                   \t');
                    fprintf(exp,'%0.3f',meanDTM);
                    
                    fprintf(exp,'\nOrder of the Taylor series                       \t'); 
                    fprintf(exp,'%0.0f',TR);
                end                
                
                fprintf(exp,'\nComputation time (s)                             \t'); 
                fprintf(exp,'%0.0f',cas);
                fprintf(exp,'\nComputation of fully normalized ALFs             \t');
                    
                if volbaALFs==1
                    fprintf(exp,'Standard forward column method');
                elseif volbaALFs==2
                    fprintf(exp,'Modified forward column method');
                elseif volbaALFs==3
                    fprintf(exp,'Extended-range arithmetic');
                end

                fprintf(exp,'\n\nExported data file contains the following columns:');
                    
                if coord==1 %Spherical coordinates
                    fprintf(exp,'\nSpherical Latitude (deg)\t');
                elseif coord==0 %Ellipsoidal coordinates
                    fprintf(exp,'\nEllipsoidal Latitude (deg)\t');
                end
                fprintf(exp,'Longitude (deg)\t');                               
                               
                for i=1:pocetpar
                    if volbapar(i)~=1
                        Par=cellstr(get(findobj('tag',sprintf('P%0.0d',i)),'string'));
                        Par=Par(volbapar(i));

                        fprintf(exp,'%s\t',char(Par));
                        
                        if volbapar(i)==2 || volbapar(i)==3 
                            fprintf(exp,'(arcsec)\t');
                        elseif volbapar(i)==4
                            fprintf(exp,'(arcsec)\tAzimuth\t(deg)\t');
                        elseif volbapar(i)==5 || volbapar(i)==11 || volbapar(i)==18
                            fprintf(exp,'(m^2.s^-2)\t');
                        elseif volbapar(i)==10 || volbapar(i)==22 || volbapar(i)==23
                            fprintf(exp,'(m)\t');
                        elseif volbapar(i)==16 || volbapar(i)==17 || volbapar(i)==19 || volbapar(i)==20 || volbapar(i)==21                        
                            fprintf(exp,'(mGal)\t');
                        elseif volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15 || volbapar(i)==24 || volbapar(i)==25
                            fprintf(exp,'(Eotvos)\t');
                        end
                    end
                end
                
                if get(findobj('tag','export'),'value')==0 && get(findobj('tag','datamat'),'value')==0
                    export_error=0;
                end
                
                if export_error==1
                    fprintf(exp,'\n\nNote: isGrafLab failed to create the data file due to the lack of memory.');
                    fprintf(exp,'\nHowever, isGrafLab created three *.mat files, which contain all the computed data.');
                end
                
                fclose(exp); 
                
                set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow; 
            end
            
            %% Display data
            %==============================================================
            
            if display_data==1                               
                
                if get(findobj('tag','export'),'value')==1 || get(findobj('tag','datamat'),'value')==1
                    fiGrid=reshape(fi,j_fi,i_fi);
                    lambdaGrid=reshape(lambda,j_fi,i_fi);
                elseif get(findobj('tag','report'),'value')==1
                    [lambdaGrid,fiGrid]=meshgrid(lambda,fi);
                else
                    fi=180/pi*(fi);
                    lambda=180/pi*(lambda);
                    [lambdaGrid,fiGrid]=meshgrid(lambda,fi);
                    [j_fi,i_fi]=size(fiGrid);  
                end
                
                clear fi lambda
                
                volbaformat=get(findobj('tag','nmin'),'userdata');
                colmap=get(findobj('tag','nmax'),'userdata');

                %Colormap
                if colmap==1
                    colmap='jet';
                elseif colmap==2
                    colmap='hsv';
                elseif colmap==3
                    colmap='hot';
                elseif colmap==4
                    colmap='cool';
                elseif colmap==5
                    colmap='spring';
                elseif colmap==6
                    colmap='summer';
                elseif colmap==7
                    colmap='autumn';
                elseif colmap==8
                    colmap='winter';
                elseif colmap==9
                    colmap='gray';
                elseif colmap==10
                    colmap='bone';
                elseif colmap==11
                    colmap='copper';
                elseif colmap==12
                    colmap='pink';
                elseif colmap==13
                    colmap='lines';
                end
                
                %Graphic format file
                if volbaformat==1
                    format='bmp';
                    dformat='-dbmp16m';
                elseif volbaformat==2
                    format='emf';
                    dformat='-dmeta';
                elseif volbaformat==3
                    format='eps';
                    dformat='-depsc';
                elseif volbaformat==4
                    format='jpeg';
                    dformat='-djpeg';
                elseif volbaformat==5
                    format='pdf';
                    dformat='-dpdf';
                elseif volbaformat==6
                    format='png';
                    dformat='-dpng';
                elseif volbaformat==7
                    format='tiff';
                    dformat='-dtiff';
                end
                
                coast=load('coast.mat'); %Loading continents  
                
                zobr=1;
      
                for i=1:pocetpar
                                       
                    if i==1
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 1st functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    elseif i==2
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 2nd functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    elseif i==3
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 3rd functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    elseif i==4
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 4th functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    end                    
                                                                              
                    if volbapar(i)~=1                                     
                        if volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15                   
                            for j=0:2
                                Pdisp=reshape(P(:,zobr+j),j_fi,i_fi);

                                if volbapar(i)==6
                                    nazov=['Disturbing tensor - Trr';...
                                        'Disturbing tensor - Tpp';...
                                        'Disturbing tensor - Tll'];
                                elseif volbapar(i)==7
                                    nazov=['Disturbing tensor - Trp';...
                                        'Disturbing tensor - Trl';...
                                        'Disturbing tensor - Tpl'];
                                elseif volbapar(i)==8
                                    nazov=['Disturbing tensor - Txx';...
                                        'Disturbing tensor - Tyy';...
                                        'Disturbing tensor - Tzz'];
                                elseif volbapar(i)==9
                                    nazov=['Disturbing tensor - Txy';...
                                        'Disturbing tensor - Txz';...
                                        'Disturbing tensor - Tyz'];
                                elseif volbapar(i)==12
                                    nazov=['Gravitational tensor - Vrr';...
                                        'Gravitational tensor - Vpp';...
                                        'Gravitational tensor - Vll'];
                                elseif volbapar(i)==13
                                    nazov=['Gravitational tensor - Vrp';...
                                        'Gravitational tensor - Vrl';...
                                        'Gravitational tensor - Vpl'];
                                elseif volbapar(i)==14
                                    nazov=['Gravitational tensor - Vxx';...
                                        'Gravitational tensor - Vyy';...
                                        'Gravitational tensor - Vzz'];
                                elseif volbapar(i)==15
                                    nazov=['Gravitational tensor - Vxy';...
                                        'Gravitational tensor - Vxz';...
                                        'Gravitational tensor - Vyz'];
                                end

                                Okno=figure('numbertitle','off','name',...
                                    char(nazov(j+1,:)),'visible','off');
                                worldmap([min(min(fiGrid)) max(max(fiGrid))],...
                                    [min(min(lambdaGrid)) max(max(lambdaGrid))]);
                                contourfm(fiGrid,lambdaGrid,Pdisp,ncolor,...
                                    'linestyle','none');
                                colormap(sprintf('%s(%d)',colmap,ncolor));
                                labelcolbar=colorbar('location','southoutside');
                                caxis([min(min(Pdisp)) max(max(Pdisp))]);

                                %Units in colorbar
                                set(get(labelcolbar,'xlabel'),'string','Eotvos');
                                
                                geoshow(coast.lat,coast.long,'color','black');
                                
                                print(Okno,sprintf('%s',dformat),...
                                    sprintf('-r%i',DPI),...
                                    sprintf('%s\\%s_%s.%s',outadresar,...
                                    outname,char(nazov(j+1,:)),format));
                            end                                                       
                            
                            zobr=zobr+3;                          
                        else                                
                            Pdisp=reshape(P(:,zobr),j_fi,i_fi);
                            
                            %Deleting azimuth if Theta has been computed
                            if volbapar(i)==4
                                P(:,zobr+1)=[];
                            end
                            
                            Nazov_okna=cellstr(get(findobj('tag','P1'),'string'));
                            Nazov_okna=Nazov_okna(volbapar(i));                            
                            Okno=figure('numbertitle','off','name',...
                                char(Nazov_okna),'visible','off');
                            worldmap([min(min(fiGrid)) max(max(fiGrid))],...
                                [min(min(lambdaGrid)) max(max(lambdaGrid))]);
                            contourfm(fiGrid,lambdaGrid,Pdisp,ncolor,...
                                'linestyle','none');
                            colormap(sprintf('%s(%d)',colmap,ncolor));
                            labelcolbar=colorbar('location','southoutside');
                            caxis([min(min(Pdisp)) max(max(Pdisp))]);
                          
                            %Units in colorbar
                            if volbapar(i)==2 || volbapar(i)==3 || volbapar(i)==4
                                set(get(labelcolbar,'xlabel'),'string','arcsec');
                            elseif volbapar(i)==5 || volbapar(i)==11 || volbapar(i)==18
                                set(get(labelcolbar,'xlabel'),'string','m^2 \cdot s^{-2}');
                            elseif volbapar(i)==10 || volbapar(i)==22 || volbapar(i)==23
                                set(get(labelcolbar,'xlabel'),'string','m');
                            elseif volbapar(i)==16 || volbapar(i)==17 || volbapar(i)==19 || volbapar(i)==20 || volbapar(i)==21                        
                                set(get(labelcolbar,'xlabel'),'string','mGal');
                            elseif volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15 || volbapar(i)==24 || volbapar(i)==25
                                set(get(labelcolbar,'xlabel'),'string','Eotvos');
                            end

                            geoshow(coast.lat,coast.long,'color','black');   
                            
                            print(Okno,sprintf('%s',dformat),...
                                sprintf('-r%i',DPI),sprintf('%s\\%s_%s.%s',...
                                outadresar,outname,char(Nazov_okna),format));
                            
                            zobr=zobr+1;
                        end  
                    end                     
                end
    
                set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;
            end 
            
            set(findobj('tag','R_text'),'userdata','');
            set(findobj('tag','ell_text'),'userdata','');
          
            clear all %#ok<CLALL>

            set(findobj('tag','hlasky'),'string',...
                    'Computation has been finished','fontsize',8,...
                    'foregroundcolor','k'); drawnow; 
                
        case('Close')
            close all
    end
end