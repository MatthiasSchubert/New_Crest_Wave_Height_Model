function [height, ParaGev, uep] = sjrndn(norm,char,Type,Para1,Para2,Para3,d0,Nlftms,NSmpls)

% This function generates random 1 hour max normalized wave heights from the Schubert-Jonathan (SJ) distribution function
% for wave and crest heigths based on the results of the AWARE project and the basin tests by DHI.
%
% The model assumes that the 1h-max wave heights result from a generalized
% extreme value distribution with a negative shape parameter, so that the
% upper end point exists. The output is not sorted.

% Theoretical details of the model can be found in:
% Schubert M., et. al (2020), On the distribution of maximum crest and wave height at intermediate water depths, Journal of Ocean Engineering

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LICENSE OF THIS FUNCTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GNU General Public License v3.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright by 
% Matrisk GmbH (https://www.matrisk.com)
% Total E&P Denmark A/S (https://dk.total.com/), 
% Shell (https://www.shell.com/)

% Written by Matthias Schubert, Matrisk GmbH & Philip Jonathan, Shell
%
% Rev. A 04.12.2017 // Error handling wrt to input. TBD before final release.
%
% Rev. B 14.12.2017 // Wave height and crest height introduced, one function for the wave and crest and normalized and non-normalized space. 
%
% Rev. C 18.12.2017 // NEW Wave height and crest height model - normalization of data by average in each WG instead of average achieved Hm0 in basin. 
%
% Rev. D 13.09.2018 // Update with revised spreading values from DHI.

%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
%%%%%%%%%%%%%%%%%%%%%
%
%   Sea-state parameters for normalized analysis:
%   
%   d0~depth (model developed for d=45m) (DEFAULT=45m) -  scalar value
%   Years - Vector of number of years considered & more than 1 event (storm) per year.

%   norm: 'real' --> gives results in the real space; unit: [m]
%           Para 1: Hs ~ Significant wave height defined as four times the standard deviation of the surface elevation in [m] (e.g. Hs=14.5)  
%           Para 2: Tp ~ Spectral peak period  in [s] (e.g. Tp=12)
%           Para 3: Sp ~ directional spreading in degree (e.g. if Sp_deg=20 then
%           Sp_dir==180*2^0.5/(pi*(20+1)^0.5)=17.68

%         'norm' --> gives results in the normalized space (normalized by Hm0); unit: [-]
%           Para 1: U ~ Ursell number (e.g. U=0.25) --- Vector same size as years (lifetimes)
%           Para 2: Stp ~ Steepness (e.g. Stp=0.0356) --- Vector same size as years (lifetimes)
%           Para 3: Sp ~ directional spreading in degree (e.g. if Sp_deg=20 then
%           Sp_dir==180*2^0.5/(pi*(20+1)^0.5)=17.68.

%   Char:  'wave' --> random numbers of wave height
%          'crest' --> random numbers of crest height

%	Type: 'Hm0MEAN'--> E[Hm0]:= the variability of Hm0 is considered in the response surface. The surface considers the variability of Hm0 in the DHI tank.
%								It is a mixture of natural variability in Hm0 coming from the spectrum used for generating the waves 
%								and possible and unquantified variability of the measurements and caused by the tank experience. 
%								These are inherent uncertainties which might not occur in reality.
%
%		  'Hm0VAR'--> Hm0:= the variability of the hourly Hm0 is not considered in the response surface, i.e. the response if fitted to the data normalized by an hourly Hm0.
%                   			For real applications it is recommended to use this option. The natural variability of Hm0 in hour is determined by the spectrum. 
%								For the JONSWAP spectrum a normal distribution with a coefficient of variation is in the range of 0.03-0.04 can be assumed. 
%								For other spectra the hourly variability might be different. The result in this case can be multiplied with (randn(Nlftms,NSmpls).*0.04+1):
%								height_new=bsxfun(@times,height,(randn(Nlftms,NSmpls).*0.04+1));
%								The natural variability due to the spectrum of Hm0 will be higher for 0.5h max and converge to zero in case of a 100h maximum.
%
%
%   Nlftms - Number of Lifetimes (epistemic uncertainty) optional, default=1
%   
%
%   NSmpls - Number of samples per lifetime (aleatoric uncertainty) optional, default=1
%	
%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT:
%%%%%%%%%%%%%%%%%%%%%
%   
%   HEIGHT 
%   All result are "Number of lifetimes (NmSm_Lftm)*realizations per lifetime (NmSm_smpls)" realizations of the normalized crest/wave
%   height.
%   In each column a number of realizations of the 1h max crest/wave height is generated. The number of realizations are the realizations within one lifetime. 
%   In each row "one lifetime" is simulated using one realization of
%   beta values. The number of columns represent the number of hours in one
%   lifetime. In each row a new set of beta is used and a "new lifetime" is
%   simulated.

%   PARAGEV
%   The output is a Number of lifetimes*3 matrix of realizations of the GEV
%   parameters. This might be of use in the case the distributions is
%   needed in an analysis.

%   UEP
%   The realization of the upper end point of the distribution is provided.
%   The variability of this value might be quite large and this value gets
%   numerically instable in the case the scale parameter is close to zero.
%   The mean value of the realizations might be stable.

%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT REMARKS:
%%%%%%%%%%%%%%%%%%%%%
%
% The model is constructed using mean value estimates for Steepness and Ursell-Number and Tp.
% In case Type 2 is chosen, the hourly variability of Hm0 conditional on the sea state parameters need to be considered outside this function.
% In Type 1 included the hourly variability of Hm0 in response surface.

% No MATLAB Toolboxes needed to be called for using this function.
% This code used a Cholesky decomposition for generating multivariate
% random numbers.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples
%           [cr,para,uep]=sjrndn('real','crest','Hm0MEAN',14.5,12,17.68,45,1,10) (only one lifetime, parameter & uep are always constant for one lifetime)
%       
%           [cr,para,uep]=sjrndn('real','crest','Hm0MEAN',14.5,12,17.68,45,10,1) (10 lifetimes)
%
%           [h,~,~]=sjrndn('real','wave','Hm0MEAN',14.5,12,17.68,45,1,10)
%
%           [cr,~,~]=sjrndn('norm','crest','Hm0MEAN',0.25,0.0356,17.68,45,1,10) 
%
%           [cr,~,~]=sjrndn('norm','crest','Hm0MEAN',0.25,0.0356,17.68,45,10,10) 
%
%           [h,~,~]=sjrndn('norm','wave','Hm0MEAN',0.25,0.0356,17.68,45,1,10) 
%
%           [h,~,~]=sjrndn('norm','wave','Hm0MEAN',0.25,0.0356,17.68,45,10,10) 
%           
%           % In the next examples Type 2 regression is chosen - variability of Hm0 has to be added
%           ouside this function, e.g. with N~(Hm0,Hm0*0.04).
%
%           [cr,~,~]=sjrndn('norm','crest','Hm0VAR',0.25,0.0356,17.68,45,1,10)
%        
%           [cr,~,~]=sjrndn('norm','crest','Hm0VAR',0.25,0.0356,17.68,45,1,10)
%
%           [h,~,~]=sjrndn('norm','wave','Hm0VAR',0.25,0.0356,17.68,45,1,10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error handling of the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~exist('d0','var')==1
            %default water depth is 45 m.
            d0 = 45; 
            warning('default water depth of 45m is used since input is missing') 
        end
        if ~exist('NSmpls', 'var')==1
            NSmpls = 1;
        end
        if ~exist('Nlftms', 'var')==1
            Nlftms = 1;
        end
        
        if nargin > 9 
        error('Too many Inputs for function crx_1hmax_norm');
        end
        
        if Para3 > 90 
        error('directional spreading should not be larger than 90[°]');
        end
        
        if d0<5 || d0>100
            warning('The input water depth is lower than 5m or more then 100m. please check your input') 
        elseif Para1 <=0 || Para2 <=0 || Para3 <=0
           error('Check input: parameters should be positive.'); 
        end
        
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
%  INPUT DATA - This data is the output of the AWARE project.
%               The model is based on the model. 
%               The response surface model is hard coded for an easy use. 
%               THIS DATA MUST NOT BE CHANGED!

if strcmp(norm,'real')
    
    Hm0=Para1;
    Tp=Para2;
    Sp=Para3;
    
    % Calculation of the local wave length using the dispersion relation (linear theory)
    Lguess=9.81*Tp^2/(2*pi);
    Lp = fzero(@(z) z - (9.81*Tp^2/(2*pi))*tanh(2*pi*d0/z),Lguess);
    %Calculation of Stp~Steepness (e.g. Stp=0.0356)
    Stp=Hm0/Lp;
    %Calculation of the wave number
    kp=2*pi/Lp;
    %Calculation of U~Ursell number (e.g. U=0.25) number
    U=Hm0/(kp^2*d0^3);
    
elseif strcmp(norm,'norm')
    
    if Para2>1
    warning('Check input - for normalized anaylsis the input is Ursell Number,Stepness and Sprading');
    end

    U=Para1;
    Stp=Para2;
    Sp=Para3;
    
else  
    error('Choose either "norm" for normalized random numbers or "real" for real space random numbers ');
end %real or normalized analysis


if strcmp(char,'wave')
    
        if strcmp(Type,'Hm0MEAN') % Type 1 E[Hm0] considered in the normalization (DEFAULT)

                %WAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Mean of the beta values of the response surface for..
                %..location parameter of the GEV (variable 'Pmodel' in regression analysis)	
                Lcbtpsi_mu=[1.48979807426707,-0.101231733049332,8.96403805822678,-0.00206401659181570,0.264540571430569,-137.429328105851,-4.69351219935953e-05,-5.16397591320618,-0.00204200820012409,0.115373969042408];
                
                Lcbtsgmh_mu=[-1.45101335586540,-0.722917488658706,-22.7647498249478,-0.0525020161452659,0.327276162990079,-56.0062671095275,0.000245562065740186,-18.4269037329312,0.0414162023874748,1.02598470056914];
                
                %..shape parameter of the GEV	
                Lcbtxi_mu=[-3.57313133424849,-2.65791497265269,97.8710855293671,-0.0386021079892585,1.44523961199852,-707.406560616691,0.00136726300545854,137.793847062836,-0.111410233969954,-0.840015253315663];
                
                %..standard error
                eps_mu=[-0.00895435831906089,0.0213984479687506,-0.379098252561749];
                
                %WAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Standard deviation of the beta values of the response surface
                %for..

                %..location parameter of the GEV     
                Lcbtpsi_std=[0.0568690747011864,0.0876230324044125,1.68268414122940,0.00259260663393709,0.0589467733991264,15.9937606306026,4.32141596527773e-05,1.30333102381380,0.00228695643974895,0.0309812247147897];
                
                %..scale parameter of the GEV	
                Lcbtsgmh_std=[0.332760537937054,0.526368044059455,10.3553661169937,0.0144441379149088,0.361194900984819,100.416525074352,0.000249828849794585,8.20435799968735,0.0128675640273077,0.180110367358209];
                
                %..shape parameter of the GEV	
                Lcbtxi_std=[2.08664785655904,3.09977739688718,65.2173924021553,0.0897262104210650,2.40969113042495,636.506356183900,0.00158884028984410,49.4203085164929,0.0747751123551964,1.22429685757956];
                
                %..standard error
                eps_std=[0.00430930030030821,0.0296039789620390,0.178825667289663];
                
                %WAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Choleski decomposition of correlation matrix xi, psi,sigma and eps

                chol_psi=[1,-0.421068677521132,-0.907278168598446,-0.808027529413277,0.386580611224819,0.742024927653311,0.633176757826785,0.336833869262241,0.168488063638170,0.644321455693038;0,0.907028758535586,-0.0570673490333067,-0.226390353599391,-0.413000321385029,0.0807268007688659,0.247281133967060,-0.726866965402535,-0.528713394297754,0.269452638297558;0,0,0.416640903487599,-0.520945003010108,-0.0680246983995070,-0.586727376427709,0.463744931484915,-0.237649211021110,0.295255285140067,0.246905656166120;0,0,0,0.156381659238604,-0.0822602109372785,-0.0503392691143572,-0.241642718415392,0.00371093635728201,-0.0190090716553566,-0.0522251583618035;0,0,0,0,0.817674790894963,0.0462066731689684,0.188732130951026,-0.0419773115947572,-0.473532962941801,-0.106292781710973;0,0,0,0,0,0.306535605547983,0.381684435337684,-0.0890066418398572,0.129635019229126,-0.625132747154144;0,0,0,0,0,0,0.288419562197914,0.449833700693538,-0.511067866942957,-0.206138665879737;0,0,0,0,0,0,0,0.299469205957408,-0.288346186827529,-0.0129270301705681;0,0,0,0,0,0,0,0,0.138435174759493,-0.0500448701079197;0,0,0,0,0,0,0,0,0,0.0360751094687323];
                
                chol_sgmh=[1.00000000000000,-0.393037161882795,-0.912282989208396,-0.785702635701719,0.443130047621356,0.764991500439664,0.592984670877045,0.280729609793727,0.124491938321350,0.623700958381946;0,0.919522587748183,-0.0763611749849878,-0.212594848859818,-0.429397226122786,0.110863569201633,0.208680218726177,-0.743022649401259,-0.524091496405614,0.273457750972431;0,0,0.402378824686271,-0.557667598858173,0.0524596571430261,-0.555988723393563,0.521469315539358,-0.185340926666902,0.183064888010099,0.236630792270775;0,0,0,0.162731827475150,-0.0685758471123332,-0.0521494064817337,-0.238487688780255,0.00813055268183245,-0.0239902269626407,-0.0487739837545585;0,0,0,0,0.782175888568813,0.0473380082721725,0.171724188586017,-0.0898380759435358,-0.422595498909450,-0.102000042177009;0,0,0,0,0,0.297343849969307,0.388733548220935,-0.0591599639031989,0.104998168995701,-0.643990364528570;0,0,0,0,0,0,0.308888577338882,0.482405845709825,-0.602940851231034,-0.220063129331543;0,0,0,0,0,0,0,0.300673714321457,-0.315759389285030,-0.0131064423149337;0,0,0,0,0,0,0,0,0.151288939716730,-0.0520214709315984;0,0,0,0,0,0,0,0,0,0.0375705767266583];
                
                chol_xi=[1.00000000000000,-0.479282007362143,-0.924917326423087,-0.752119680889606,0.444345149289661,0.771067014046052,0.548170682774020,0.393413711910104,0.0942503003363558,0.557806940408026;0,0.877660958126152,-0.0463337605357872,-0.215187897598455,-0.423424151942497,0.0646158409842496,0.144688429964785,-0.756818117028662,-0.360714467088603,0.269033612973033;0,0,0.377334231043226,-0.598339484457626,-0.0182418184428823,-0.552447597855697,0.488708887040520,-0.121349566485107,0.194502915726073,0.332208684406704;0,0,0,0.173205126034588,-0.121694917255820,-0.0614584960184551,-0.260191391905801,0.0289435448540866,-0.00366612797500214,-0.0390507612074105;0,0,0,0,0.779824954094522,0.0513361378549515,0.196295369950519,-0.0283863721496432,-0.568884407116643,-0.115297780340555;0,0,0,0,0,0.299448758195394,0.479247030695928,-0.0317576393635068,0.0798643335607125,-0.669277532169366;0,0,0,0,0,0,0.322224510849584,0.423905034206355,-0.598290124964782,-0.199322448135306;0,0,0,0,0,0,0,0.274551873306681,-0.330996665585527,-0.0148847286334071;0,0,0,0,0,0,0,0,0.160123335915027,-0.0460432114624301;0,0,0,0,0,0,0,0,0,0.0358858796758734];
                
                chol_eps=[1,0.285761515922544,0.343268710932846;0,0.958300764905074,0.466248809904807;0,0,0.815339585299802];
                

        elseif strcmp(Type,'Hm0VAR') % Type 2 Hm0 in each hour is considered in the normalization.

                %WAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Mean of the beta values of the response surface for..
                %..location parameter of the GEV (variable 'Pmodel' in regression analysis)	
                Lcbtpsi_mu=[1.50257113603886,-0.0605403819717037,8.67542639615047,-0.00249741338374870,0.241622921663811,-134.786681354542,-4.37601095910506e-05,-5.98174523829633,-0.00189531595738265,0.122401949193713];
                
                Lcbtsgmh_mu=[-1.66026155969105,-0.552540575204619,-20.6080205777377,-0.0455669279826572,-0.0446025241470169,-28.6356389133861,0.000263780578292916,-12.6923185217040,0.0350787418817457,0.804292521905260];
                
                %..shape parameter of the GEV	
                Lcbtxi_mu=[-1.51162834843763,1.58654768698325,4.11375430468567,-0.0580681661256932,-1.79849417301255,426.091460398549,0.00146174387998523,90.6613394326347,-0.109318630145308,-0.789606418120764];
                
                %..standard error
                eps_mu=[-0.00627506773615667,0.0163985996077762,-0.249858067469834];
                
                %WAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Standard deviation of the beta values of the response surface
                %for..

                %..location parameter of the GEV     
                Lcbtpsi_std=[0.0544321668508869,0.0861321991583864,1.61946423831097,0.00242988472481378,0.0531746229423002,15.5522269253287,3.97504370536150e-05,1.24678758259374,0.00209181918611976,0.0288822401651744];
                
                %..scale parameter of the GEV	
                Lcbtsgmh_std=[0.335588046720347,0.561912459921343,10.1459610603365,0.0149026375023901,0.355717287061905,97.9483009430339,0.000253570957054564,8.35384019920551,0.0137016767037427,0.192058983767931];
                
                %..shape parameter of the GEV	
                Lcbtxi_std=[1.83305935561377,2.61704353832367,59.0056073722300,0.0791928363598361,1.97315962292515,584.535992388722,0.00137917957551661,41.4649992995236,0.0766225336721760,1.05545652335530];
                
                %..standard error
                eps_std=[0.00411227342642945,0.0277614041636772,0.132808834497159];
                
                %WAVE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Choleski decomposition of correlation matrix xi, psi,sigma and eps

                chol_psi=[1.00000000000000,-0.363252340659423,-0.916622789760271,-0.819046639566472,0.377936456680100,0.763009079040738,0.649941246141321,0.256199140297723,0.204153256661347,0.651820941318064;0,0.931690794741179,-0.0696206288958457,-0.226405048180172,-0.468670932710873,0.101309257969289,0.260942860372706,-0.776497040929662,-0.566749651130493,0.269538586445854;0,0,0.393644038852674,-0.505198253989244,0.00425872979346523,-0.564933803411415,0.453804368490344,-0.197463539733935,0.222194761411211,0.265105699462337;0,0,0,0.150592431880309,-0.0483337498636432,-0.0449754225115193,-0.221425997811705,0.00765168397150690,-0.0484660596271628,-0.0597205958325673;0,0,0,0,0.796967567338016,0.0477117038835095,0.181222616903443,-0.0642161442856411,-0.424402794434545,-0.106727179830282;0,0,0,0,0,0.290007210249459,0.374563305754769,-0.0905213015614940,0.134725697590647,-0.612014555812685;0,0,0,0,0,0,0.285268322674769,0.433050823465124,-0.525265592530836,-0.197198424604112;0,0,0,0,0,0,0,0.304159683750507,-0.304447327959961,-0.0142377757653128;0,0,0,0,0,0,0,0,0.136146781415477,-0.0479092462690175;0,0,0,0,0,0,0,0,0,0.0359596048048081];
                
                chol_sgmh=[1.00000000000000,-0.396507299861629,-0.913171994014409,-0.797285401136377,0.393229756251549,0.752748446679845,0.603028900731529,0.321890211493941,0.166034150031324,0.601436497446408;0,0.918031568714519,-0.0621683594977586,-0.240312319068604,-0.427381665797298,0.0656465051964536,0.229753957113421,-0.748464948465376,-0.585417764442395,0.317264432305286;0,0,0.402805169312792,-0.531301461049086,-0.00423481096260548,-0.567933634852341,0.486424359884041,-0.169482352532872,0.199877442968943,0.232264138964210;0,0,0,0.155899762430169,-0.0419712583340681,-0.0492137077925745,-0.227160555037526,-0.00888187666269762,-0.0278372480795253,-0.0471950333191079;0,0,0,0,0.812979551027356,0.0112586704862235,0.129401129759051,-0.0810770372976550,-0.428755168819448,-0.0220237518092125;0,0,0,0,0,0.322432864460449,0.431907794935334,-0.0138241110257791,0.0475217814867617,-0.658189579048386;0,0,0,0,0,0,0.303429881044817,0.456291001781469,-0.541007238980836,-0.210191349512101;0,0,0,0,0,0,0,0.304002774297584,-0.302792339789166,-0.0133599488202803;0,0,0,0,0,0,0,0,0.136130266070655,-0.0467814606890593;0,0,0,0,0,0,0,0,0,0.0346042453548803];
                
                chol_xi=[1.00000000000000,-0.442808414009568,-0.920477953895830,-0.752180730407967,0.397771791581403,0.783618456461941,0.556060373567853,0.322172237918877,0.153924717409932,0.568504440206332;0,0.896616254861762,-0.0721310173434410,-0.174433080250567,-0.334481004595002,0.0771673230239205,0.133256453511064,-0.721933320826364,-0.436305502428156,0.278123504856592;0,0,0.384080008238840,-0.611698981313152,-0.0233996885264824,-0.538533392944262,0.525445261892256,-0.184258985612538,0.226550323229038,0.319478861271574;0,0,0,0.172109283821903,-0.0506036990548885,-0.0141622954164136,-0.202767722999004,0.00607805505793309,-0.0304203457797182,-0.126139473053338;0,0,0,0,0.852520838223950,0.0745391414667546,0.266855014917216,-0.00252485873957150,-0.556767686585268,-0.182824355590990;0,0,0,0,0,0.290193813530552,0.404542092371980,-0.136518240294573,0.178943740313920,-0.628548638963308;0,0,0,0,0,0,0.347801958480463,0.475485105669737,-0.536374947039631,-0.220784825384048;0,0,0,0,0,0,0,0.310321321059594,-0.294531925373332,-0.0173911869758309;0,0,0,0,0,0,0,0,0.131285892729234,-0.0512911651407837;0,0,0,0,0,0,0,0,0,0.0359832653428755];
                
                chol_eps=[1,0.214324261501873,0.358571812176679;0,0.976762566303438,0.546664865316172;0,0,0.756692659235693];
                
        else %choose type
             error('Choose either "Hm0MEAN" or "Hm0VAR" for the option Type. Hm0MEAN should be taken as default.');
        end 
 
elseif strcmp(char,'crest')
   
        if strcmp(Type,'Hm0MEAN') % Type 1 E[Hm0] considered in the normalization (DEFAULT)
 

            %CREST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Mean of the beta values of the response surface for..
            %..location parameter of the GEV 	
            Lcbtpsi_mu=[0.835602772827991,0.196930349321298,8.88694345216497,-0.00139849473093344,0.00520503428185160,-106.325572720999,-1.33382935933788e-05,0.476106762609218,-0.00564752436363829,0.0267230220769566];
            
	    %.scale parameter of the GEV
            Lcbtsgmh_mu=[-2.59398920521512,0.444376697890503,-3.08646521022753,-0.0225739665238025,-0.184761314071880,-198.623965715679,-0.000176972592312079,-15.9543108143504,0.0186614276234371,0.724028022372064];
            
            %..shape parameter of the GEV	
            Lcbtxi_mu=[-7.02101235370273,-2.95459587413733,189.418740277125,0.116923383013318,5.31120642603432,-1212.01475717851,-0.000364945008739898,85.4748767500036,-0.111141481498853,-2.87978983344630];
            
            %..standard error
            eps_mu=[0.00141582583627388,0.00302388840544803,-0.112438233319659];
            
            %CREST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Standard deviation of the beta values of the response surface
            %for..

            %..location parameter of the GEV     
            Lcbtpsi_std=[0.0402067013429670,0.0686433691801948,1.15108205219201,0.00178401548934329,0.0463561739788465,10.7079029392985,2.94990354356842e-05,1.01254397395073,0.00161891624048228,0.0205263930511906];
            
            %..scale parameter of the GEV	
            Lcbtsgmh_std=[0.281243578713025,0.508867135989823,8.34187861149519,0.0122529990816789,0.323122426506744,80.8031858210346,0.000206590151433163,7.56376757466452,0.0124122320488671,0.152371625760973];
            
            %..shape parameter of the GEV	
            Lcbtxi_std=[2.52577427759224,3.52891841328450,76.6597920922429,0.109386310213763,2.35369887020221,725.659447493706,0.00167380439009607,51.8254120991350,0.0992231991017831,1.47840798062379];
            
            %..standard error
            eps_std=[0.00371654680853546,0.0277151521034850,0.174790246030054];
            
            %CREST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Choleski decomposition of correlation matrix xi, psi,sigma and eps

            chol_psi=[1.00000000000000,-0.469972621573463,-0.921311415236750,-0.819492288926081,0.491352027195746,0.767230032449070,0.662261599733619,0.348030907741430,0.183814350485866,0.658616032756331;0,0.882680992755234,-0.0311689257585634,-0.263440356158488,-0.358556474387276,0.0502513214062130,0.277701373063293,-0.760496979298022,-0.561996064645554,0.310848745934598;0,0,0.387574217694517,-0.487467924853733,-0.0355097661613202,-0.572341967701355,0.451687385891788,-0.148151495340176,0.202619175353425,0.246561821840983;0,0,0,0.146309908675314,-0.0283629839219857,-0.0278615689542653,-0.195873168741414,-0.00192412347317713,-0.0345507448389050,-0.0969337390138021;0,0,0,0,0.792429831402004,0.0646106447653574,0.191611370759045,-0.0751327914724159,-0.430507005228876,-0.135826691089137;0,0,0,0,0,0.276236767232898,0.344161970911773,-0.0840206012647582,0.123384446186039,-0.575862765088166;0,0,0,0,0,0,0.294518842324960,0.446116246660171,-0.561285958120800,-0.213438820135583;0,0,0,0,0,0,0,0.258538508905545,-0.264753166391590,-0.0139320405299247;0,0,0,0,0,0,0,0,0.149760657980719,-0.0476338649128197;0,0,0,0,0,0,0,0,0,0.0363739713327748];
            
            chol_sgmh=[1.00000000000000,-0.516739034823848,-0.913788003830410,-0.787771334788579,0.490615421818487,0.738863049723602,0.586211739402403,0.388668720128301,0.242536808635776,0.613320243739811;0,0.856142961128175,-0.0435300015171836,-0.265576787457950,-0.337092588614602,0.0772799578957644,0.311741785442860,-0.718376603433627,-0.567503555162982,0.266995529496083;0,0,0.403852229192250,-0.530047071294117,0.00905538858294409,-0.595893711031644,0.453851566904747,-0.176439400732224,0.198855863591319,0.294351746346590;0,0,0,0.167138853237260,-0.0559511383475858,-0.0457165549157217,-0.238367349213348,5.07016203185426e-05,-0.0205508460847083,-0.0713906983283646;0,0,0,0,0.801531387177195,0.0520029648265827,0.194103215290704,-0.0721853747252016,-0.403332047759752,-0.117182564680843;0,0,0,0,0,0.297027903183497,0.395124935625706,-0.0485358415766373,0.0764672121387420,-0.626413800513953;0,0,0,0,0,0,0.320269956292200,0.458825678917646,-0.552227623804600,-0.224328037170842;0,0,0,0,0,0,0,0.289228900795865,-0.292125082227134,-0.0121794892002337;0,0,0,0,0,0,0,0,0.142594642722898,-0.0526014480717153;0,0,0,0,0,0,0,0,0,0.0380494410745839];
            
            chol_xi=[0.999999999999999,-0.514932251229490,-0.932744166605246,-0.805388748290747,0.362951677479828,0.767710920995359,0.602591624374805,0.385434800342065,0.292110912489504,0.666188733124594;0,0.857230877094222,-0.0651578346724121,-0.152687140598679,-0.333902964292396,0.103640385129644,0.177874521913652,-0.728298275829765,-0.476123045071736,0.180595727960929;0,0,0.354602278961499,-0.550375281952612,-0.126883574851318,-0.567435601246549,0.454394712082635,-0.173721414852855,0.277388847098151,0.381696532361691;0,0,0,0.158501262575559,0.0176110392447621,-0.0201039646065504,-0.188901587098758,0.0320793075778672,-0.0996392979623324,-0.120308585259980;0,0,0,0,0.860444826808203,0.0940206931095508,0.304015773987188,0.0135520595862845,-0.515198993158879,-0.214657543466377;0,0,0,0,0,0.262014103946531,0.421040486963276,-0.0994163521555472,0.130691778032086,-0.536200258962443;0,0,0,0,0,0,0.305589499516806,0.451649781391396,-0.488527414256320,-0.166046801164314;0,0,0,0,0,0,0,0.275242855910320,-0.256588230010846,-0.0149347108409218;0,0,0,0,0,0,0,0,0.118737529476262,-0.0347453373329296;0,0,0,0,0,0,0,0,0,0.0286522858829816];
            
            chol_eps=[1.00000000000000,0.314767826606544,0.312066951890850;0,0.949168696983519,0.453011660664709;0,0,0.835101582347534];
            
         elseif strcmp(Type,'Hm0VAR') % Type 2 Hm0 in each hour is considered in the normalization.
    
            %CREST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Mean of the beta values of the response surface for..
            %..location parameter of the GEV 	
            Lcbtpsi_mu=[0.844941619916275,0.213763010177831,8.59299241397393,-0.00161548780617506,0.00580974362119479,-101.163250514980,1.08167063384216e-06,0.540151058052228,-0.00639354384158727,0.0195314031704004];	

            %..scale parameter of the GEV		
            Lcbtsgmh_mu=[-2.91798701785471,0.764812631220784,2.05016623748413,-0.00975666232719819,-0.255264173192931,-216.382579983198,-0.000303322034756579,-16.1128392392643,0.00507473417477184,0.583141773546287];

            %..shape parameter of the GEV	
            Lcbtxi_mu=[-6.18994423898820,-3.77345664581781,129.528679443770,0.175843547901511,7.87153404062685,-514.307031814883,-0.00126264958144635,85.7451168254957,-0.146478937906443,-3.08520776327232];

            %..standard error
            eps_mu=[0.00251794954894934,-0.0113638003856248,-0.00147056347776003];

            %CREST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Standard deviation of the beta values of the response surface
            %for..

            %..location parameter of the GEV     
            Lcbtpsi_std=[0.0379348653276636,0.0625102171675321,1.09951536904150,0.00172080841121614,0.0419291657320275,10.5227854535279,2.86740289204109e-05,0.973189366631175,0.00147675006011275,0.0198334174575715];
            
            %..scale parameter of the GEV	
            Lcbtsgmh_std=[0.285161679292978,0.498308824876765,8.50977577551321,0.0128368563784688,0.340828208493052,82.3702061781608,0.000217327636065032,7.70990726726322,0.0119501077273551,0.157097205457748];
            
            %..shape parameter of the GEV	
            Lcbtxi_std=[2.47537042799661,3.56621770511297,75.0653488177880,0.116532415293064,2.31464687869880,712.733129075886,0.00170290721064401,49.4586972174500,0.107249024694867,1.54004199951951];
            
            %..standard error
            eps_std=[0.00342310330373490,0.0264692027722988,0.150149698865849];
            
            %CREST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Choleski decomposition of correlation matrix xi, psi,sigma and eps

            chol_psi=[1.00000000000000,-0.418717703490922,-0.920749773894801,-0.821081276907002,0.444878360685817,0.759137286281716,0.678666240669029,0.294215426389464,0.169560406423134,0.646453972166209;0,0.908116448911310,-0.0315810142967023,-0.257928812894677,-0.395356171270250,0.0496595103417033,0.259623684260977,-0.780311549823541,-0.510559047100956,0.309552922152542;0,0,0.388873364231422,-0.485548093863026,-0.0579989849211725,-0.584661782253629,0.410420992548705,-0.156670891564964,0.232173282416448,0.303135002693860;0,0,0,0.153431785290702,-0.0458867747775626,-0.0688381270222733,-0.259369867974940,-0.000434140969149538,-0.0278586719573716,-0.0144629166671654;0,0,0,0,0.800192016756389,0.0532786467972196,0.179903622208699,-0.0708885429703746,-0.440155297075913,-0.117265084036945;0,0,0,0,0,0.268025767501907,0.348901859509244,-0.0662446462474571,0.0993087115824485,-0.576173271971018;0,0,0,0,0,0,0.286691011292832,0.444076702396361,-0.586313399796450,-0.211276874178928;0,0,0,0,0,0,0,0.270901375138938,-0.289628851785155,-0.0149628962577158;0,0,0,0,0,0,0,0,0.157006398017150,-0.0480938833888627;0,0,0,0,0,0,0,0,0,0.0356746089815376];
            
            chol_sgmh=[1.00000000000000,-0.476065543513694,-0.905959775927097,-0.786186817693832,0.467543519467289,0.728673714952514,0.600271692313830,0.336919003423804,0.203537505756181,0.623642360558713;0,0.879409801104702,-0.0441062143496646,-0.258378667600892,-0.398184504574434,0.0809410725079890,0.281367141079160,-0.736076394930247,-0.511494076960242,0.269245401974232;0,0,0.421060003156161,-0.537994417203425,-0.0181520059571719,-0.607949035656272,0.482167731072310,-0.177128934724410,0.225618690032685,0.269596465106051;0,0,0,0.160351984306974,-0.0104071311560362,-0.0355791303376779,-0.211925425488697,-0.0171147130006057,-0.0368614882045247,-0.0887854795540473;0,0,0,0,0.788932414102617,0.0503026644841093,0.174505066711622,-0.0841701781984257,-0.431415173002018,-0.100605728160317;0,0,0,0,0,0.298470932241958,0.394869641988680,-0.0593579435938867,0.0920882672063272,-0.625736286454402;0,0,0,0,0,0,0.311021330457175,0.469965266039978,-0.582522167102173,-0.228141308580710;0,0,0,0,0,0,0,0.285541971330010,-0.296242414970350,-0.0118674173647240;0,0,0,0,0,0,0,0,0.151633724263848,-0.0527671324982080;0,0,0,0,0,0,0,0,0,0.0370260581597387];
            
            chol_xi=[1.00000000000000,-0.493641985368765,-0.906372125271390,-0.776862934724185,0.260575047705519,0.704095760023597,0.552367131353382,0.371385567815605,0.330308512952038,0.697312399184439;0,0.869665217357339,-0.0485464981738791,-0.145733398974139,-0.417464291930811,0.0541029932284027,0.151626862822154,-0.701064892548891,-0.490054852272645,0.197162968550280;0,0,0.419681793798682,-0.593159516611738,-0.144527382736737,-0.654565631077354,0.544948737096400,-0.178274899215614,0.261177376076556,0.388485092926195;0,0,0,0.152995244789185,0.0422782507102221,-0.0209145298144253,-0.195972781438859,0.0393832906634729,-0.109399385413153,-0.115671197311042;0,0,0,0,0.857408067743782,0.0872343155792650,0.282920272668255,0.0696580180080964,-0.530407985667666,-0.189978082296150;0,0,0,0,0,0.254595007965366,0.406767896495336,-0.0977561434506929,0.116331281327831,-0.495448569014005;0,0,0,0,0,0,0.301697509498145,0.473845551750286,-0.444421933306316,-0.164226465797928;0,0,0,0,0,0,0,0.313543750075121,-0.258581898908425,-0.0132511332441907;0,0,0,0,0,0,0,0,0.106393746104850,-0.0333337612006195;0,0,0,0,0,0,0,0,0,0.0276345424311390];

            chol_eps=[1.00000000000000,0.154520673063155,0.259136707954790;0,0.987989555408411,0.504177997110439;0,0,0.823803807845084];
            
        else
             error('Choose either "Hm0MEAN" or "Hm0VAR" for the option Type. Hm0MEAN should be taken as default.');
        end %choose type
    
else
    
error('Choose either "crest" for crest or "wave" for wave as second input for this function');

end %wave or crest
	   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random number simulation from the rgression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulate multivariate distribution of beta values and the standard error
    %[~,ia]=unique(simuyears); ia(end+1)=length(simuyears)+1;
    NmSm_Lftms=Nlftms;%length(ia)-1;
  
            % Generate iid standard normal distributed values
            r_psi=randn(NmSm_Lftms,size(Lcbtpsi_mu,2));
            r_sgmh=randn(NmSm_Lftms,size(Lcbtsgmh_mu,2));
            r_xi=randn(NmSm_Lftms,size(Lcbtxi_mu,2));
            r_eps=randn(NmSm_Lftms,size(eps_mu,2));
    
                    %multivariate standard normal distributed numbers
                    %by using the cholesky decomposition
                    r_psi=r_psi*chol_psi;
                    r_sgmh=r_sgmh*chol_sgmh;
                    r_xi=r_xi*chol_xi;
                    r_eps=r_eps*chol_eps;
                    
                    %transform to multivariate normal distributed numbers of the beta values
                    btsm_Gpsi=bsxfun(@times,r_psi,Lcbtpsi_std);
                    btsm_Gpsi=bsxfun(@plus,btsm_Gpsi,Lcbtpsi_mu);
                    
                    btsm_Gsgm=bsxfun(@times,r_sgmh,Lcbtsgmh_std);
                    btsm_Gsgm=bsxfun(@plus,btsm_Gsgm,Lcbtsgmh_mu);
                    
                    btsm_Gxi=bsxfun(@times,r_xi,Lcbtxi_std);
                    btsm_Gxi=bsxfun(@plus,btsm_Gxi,Lcbtxi_mu);
                    
                    Beps=bsxfun(@times,r_eps,eps_std);
                    Beps=bsxfun(@plus,Beps,eps_mu);
                  
 
%Calculate response surface for GEV parameters
        Gpsi=btsm_Gpsi(:,1)+ btsm_Gpsi(:,2).*U + btsm_Gpsi(:,3).*Stp + btsm_Gpsi(:,4).*Sp + btsm_Gpsi(:,5).*U.^2 + btsm_Gpsi(:,6).*Stp.^2 + btsm_Gpsi(:,7).*Sp.^2 + btsm_Gpsi(:,8).*U.*Stp + btsm_Gpsi(:,9).*U.*Sp + btsm_Gpsi(:,10).*Stp.*Sp + Beps(:,1);
        Gsgm=btsm_Gsgm(:,1)+ btsm_Gsgm(:,2).*U + btsm_Gsgm(:,3).*Stp + btsm_Gsgm(:,4).*Sp + btsm_Gsgm(:,5).*U.^2 + btsm_Gsgm(:,6).*Stp.^2 + btsm_Gsgm(:,7).*Sp.^2 + btsm_Gsgm(:,8).*U.*Stp + btsm_Gsgm(:,9).*U.*Sp + btsm_Gsgm(:,10).*Stp.*Sp + Beps(:,2);
        Gxi=  btsm_Gxi(:,1)+ btsm_Gxi(:,2).*U  + btsm_Gxi(:,3).*Stp  + btsm_Gxi(:,4).*Sp  + btsm_Gxi(:,5).*U.^2  + btsm_Gxi(:,6).*Stp.^2  + btsm_Gxi(:,7).*Sp.^2  + btsm_Gxi(:,8).*U.*Stp  + btsm_Gxi(:,9).*U.*Sp  + btsm_Gxi(:,10).*Stp.*Sp  + Beps(:,3);
   
          
%Tansform the response surface for scale and location parameter according
%to the response model developped in WS04 in the AWARE project.
         Gsgm=2.*exp(Gsgm)./(exp(Gsgm)+1); %supported only on [0 2]
         Gxi=-0.5.*exp(Gxi)./(exp(Gxi)+1); %supported only on [0 -0.5]


%Set GEV parameters for each realization in one lifetime(colums)-drawn from
%the same GEV parameter.
          Gxi=repmat(Gxi,1,NSmpls);
          Gsgm=repmat(Gsgm,1,NSmpls);
          Gpsi=repmat(Gpsi,1,NSmpls);
             
%Froude scaling factor for taking into account different legth and time
%scales when moving from the reference depth to any other depth of intrest.
         tau=(d0/45)^0.5;
         
%generate scaled uniform distributed numbers
         u = rand(NmSm_Lftms,NSmpls).^tau;
         
%Generate gev distuted 1h max values
         r = expm1(-Gxi.*log(-log(u)))./Gxi;
                 
        if strcmp(norm,'norm') % norm space
             if nargout == 1 || nargout == 0
                 %Wave height (normalized by Hm0)
                 height= Gpsi + Gsgm.*r;
             elseif nargout == 2
                 %Wave height (normalized by Hm0)
                 height= Gpsi + Gsgm.*r;
                 %Parameter of the GEV
                 ParaGev=[Gxi,Gsgm,Gpsi];
             elseif nargout == 3
                  %Wave height (normalized by Hm0)
                  height= Gpsi + Gsgm.*r;
                  %Parameter of the GEV
                  ParaGev.xi=Gxi;
				  ParaGev.sgm=Gsgm;
				  ParaGev.mu=Gpsi;
                  %Calculate samples of the upper end point for different lifetimes     
                  uep=Gpsi-Gsgm./Gxi;
             else
                error('Too many outputs for function (max 3: height, Paramters, upper end point');
             end %nargout 
             
        else % real space
             if nargout == 1 || nargout == 0
                 %Wave height (normalized by Hm0)
                 height= (Gpsi + Gsgm.*r).*Hm0;
             elseif nargout == 2
                 %Wave height (normalized by Hm0)
                 height= (Gpsi + Gsgm.*r).*Hm0;
                 %Parameter of the GEV
                 ParaGev=[Gxi,Gsgm.*Hm0,Gpsi.*Hm0];
             elseif nargout == 3
                  %Wave height (normalized by Hm0)
                  height= (Gpsi + Gsgm.*r).*Hm0;
                  %Parameter of the GEV
                  ParaGev.xi=Gxi;
				  ParaGev.sgm=Gsgm.*Hm0;
				  ParaGev.mu=Gpsi.*Hm0;
                  %Calculate samples of the upper end point for different lifetimes     
                  uep=(Gpsi-Gsgm./Gxi).*Hm0;
             else
                error('Too many outputs for function (max 3: height, Paramters, upper end point');
             end %nargout   
             
        end %(norm, real space)
        
end %function

