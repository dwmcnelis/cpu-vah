#include <math.h> // for math functions

#include "edu/osu/rhic/trunk/hydro/AnisotropicDistributionFunctions.h"

double transversePressureHat(double e, double p, double pl) {
	double a = pl/e;
	double a2 = a*a;
	double a3 = a2*a;
	double a4 = a3*a;
	double a5 = a4*a;
	double a6 = a5*a;
	double a7 = a6*a;
	double a8 = a7*a;
	double a9 = a8*a;
	double a10 = a9*a;
	double a11 = a10*a;

	PRECISION Rbar0 = (0.015267955823446243 + 7.725572805021035*a + 421.0063884634789*a2 + 3422.877939650926*a3 - 5785.670846299543*a4 - 12261.66452089229*a5 + 
     31491.409484673808*a6 - 22737.05146992673*a7 + 5441.373392185447*a8)/
   (0.05470696094814806 + 14.505878005231883*a + 522.6643024173569*a2 + 2731.7776413939037*a3 - 6161.1991042880445*a4 - 
     3989.4375208972588*a5 + 15008.260526258282*a6 - 10243.036679405379*a7 + 2116.74060159494*a8);

	return (e-pl-Rbar0*(e-3*p))/2;
}

double Rbar0_fun(double a) {
	double a2 = a*a;
	double a3 = a2*a;
	double a4 = a3*a;
	double a5 = a4*a;
	double a6 = a5*a;
	double a7 = a6*a;
	double a8 = a7*a;
	double a9 = a8*a;
	double a10 = a9*a;
	double a11 = a10*a;

	PRECISION Rbar0 = (0.015267955823446243 + 7.725572805021035*a + 421.0063884634789*a2 + 3422.877939650926*a3 - 5785.670846299543*a4 - 12261.66452089229*a5 + 
     31491.409484673808*a6 - 22737.05146992673*a7 + 5441.373392185447*a8)/
   (0.05470696094814806 + 14.505878005231883*a + 522.6643024173569*a2 + 2731.7776413939037*a3 - 6161.1991042880445*a4 - 
     3989.4375208972588*a5 + 15008.260526258282*a6 - 10243.036679405379*a7 + 2116.74060159494*a8);

	return Rbar0;
}

double Rbar0P_fun(double a) {
	double a2 = a*a;
	double a3 = a2*a;
	double a4 = a3*a;
	double a5 = a4*a;
	double a6 = a5*a;
	double a7 = a6*a;
	double a8 = a7*a;
	double a9 = a8*a;
	double a10 = a9*a;
	double a11 = a10*a;
	double a12 = a11*a;
	double a13 = a12*a;
	double a14 = a13*a;
	double a15 = a14*a;

	double Rbar0P = (2.473173363908116e-10 - 4.899839370307281e-6*a + 155055.91462124084*a10 - 275435.45350226434*a11 + 350689.68825705117*a12 - 
     299725.38986957155*a13 + 151477.08809203724*a14 - 33196.47417939176*a15 - 0.004301975027942015*a2 - 0.14858206981041563*a3 + 
     6.249255189587875*a4 - 92.79641927240235*a5 + 807.175057749925*a6 - 4760.015905266286*a7 + 20324.533122685436*a8 - 
     64758.869552496515*a9)/(0.00008222793468208523 + 0.03411917870833943*a + 4.895969276094396e6*a10 - 8.84162305829353e6*a11 + 
     1.1445063656613324e7*a12 - 9.918442713390596e6*a13 + 5.065882388219598e6*a14 - 1.1181016364928822e6*a15 + 0.23871740573818725*a2 - 
     23.50912574236691*a3 + 417.4953123877312*a4 - 4234.215775452717*a5 + 29824.022790048104*a6 - 157419.8447785501*a7 + 
     641300.6529027821*a8 - 2.0248032895288002e6*a9);

	return Rbar0P;
}

void secondOrderTransportCoefficientsZ(double e, double p, double pl, double cs2, double T,
double *beta_lPi, double *delta_lPi, double *lambda_piPi, double *beta_PiPi, double *delta_PiPi, double *lambda_Pipi) {
	double a = pl/e;
	double a2 = a*a;
	double a3 = a2*a;
	double a4 = a3*a;
	double a5 = a4*a;
	double a6 = a5*a;
	double a7 = a6*a;
	double a8 = a7*a;
	double a9 = a8*a;
	double a10 = a9*a;
	double a11 = a10*a;

//======================================================================
	// 14 moment \gamma functions
	//======================================================================
	double oneThirdMinusCs2 = 1.0/3.0-cs2;
	if(oneThirdMinusCs2==0.0) oneThirdMinusCs2=1.e-16;
	double meqHat2 = 36*oneThirdMinusCs2;
	double z2Logz = -2*(1-2/9*(e-3*p)/(oneThirdMinusCs2*(e+p))+0.02536*meqHat2);
	//---------------------------------------------------------------
	double m2g0m2000Z2LogZTerm = (-0.0007055993910846846 - 0.13159795793270815*a - 19.558366627527693*a2 - 489.8511366459512*a3 + 3349.0713852568874*a4 - 9839.72133209246*a5 + 
     15933.372928299428*a6 - 14518.978863672039*a7 + 6928.128191821795*a8 - 1342.330347385699*a9)/
   (0.00006248831061655302 + 0.015846410481313772*a + 2.5994701410163885*a2 + 95.00566850076045*a3 - 661.7123254379953*a4 + 1998.3169835920473*a5 - 
     3282.710432765608*a6 + 2999.3034910760416*a7 - 1423.5208949608364*a8 + 272.7021053579765*a9);

	double m2g0m2000Z2Term = (-0.01621048126789623 - 8.709181710590503*a - 15019.455310761608*a10 - 531.5483974451457*a2 - 354.9239710072709*a3 + 5026.475834641602*a4 - 
     2384.5421001124564*a5 - 37024.81099126616*a6 + 105944.59458460551*a7 - 125533.00421218756*a8 + 69885.83572265343*a9)/
   (0.0010706138978078513 + 0.9561088369582801*a + 1974.0394284975284*a10 + 71.58921949833918*a2 - 38.351803940866326*a3 - 413.8216821554177*a4 - 
     62.91046422360702*a5 + 5786.634119926802*a6 - 15619.500395846926*a7 + 17898.028788727293*a8 - 9596.659895556504*a9);	

	double m2g0m2000 = m2g0m2000Z2LogZTerm * z2Logz + m2g0m2000Z2Term * meqHat2;

	//---------------------------------------------------------------
	double g0m2040DividedBym2ConstantTerm = (0.0002022034619596147 - 0.3465776971894224*a - 3669.180338013278*a10 - 0.46483928190644463*a2 + 305.8182449301262*a3 - 3094.0393881108607*a4 + 
     13530.068145958274*a5 - 33008.96608710712*a6 + 47814.11155495329*a7 - 40754.13233671886*a8 + 18877.099618956072*a9)/
   (0.3231726678495768 + 6.510257837715722*a - 146.6399433638468*a10 - 66.42749605402405*a2 + 278.2072846508577*a3 - 659.8881404723605*a4 + 
     1009.6856692105162*a5 - 1129.2454150937733*a6 + 1136.3819011323233*a7 - 1029.3638914685152*a8 + 600.4581631072344*a9);

	double g0m2040DividedBym2Zm2Term = (-5.428668780471906e-6 + 0.07087180484499124*a + 14430.23620485441*a10 - 3.6056344041371626*a2 - 75.98232257237731*a3 + 314.3389560219402*a4 - 
     4860.077726215435*a5 + 28903.651438438494*a6 - 74655.13722955724*a7 + 95355.62148623986*a8 - 59409.22212071206*a9)/
   (0.06525900311126431 + 1.7911514139670908*a + 284.5874952982027*a10 - 1.0324279865217276*a2 + 81.69891535142081*a3 - 259.5916250456495*a4 - 
     59.79985417459323*a5 + 2029.823700195329*a6 - 4622.384249181224*a7 + 4479.931023062334*a8 - 1935.0883554632262*a9);

	double g0m2040DividedBym2LogZTerm = (1.5260366926367532e-7 - 0.00024895619154304205*a - 16.08135726710075*a10 + 2.635895081202281*a11 - 0.0027342593291142404*a2 + 0.20899592202076162*a3 - 
     2.1908540237701564*a4 + 10.836996467643049*a5 - 31.053876569731955*a6 + 55.16397490084364*a7 - 61.60172864783581*a8 + 42.08494163892432*a9)/
   (0.0000605632612588605 + 0.0012923594925415354*a - 0.18379144108413883*a10 + 0.038614228316574885*a11 - 0.016304038454467058*a2 + 0.07975759323101335*a3 - 
     0.21819092742589902*a4 + 0.37685396820303946*a5 - 0.45029351728017364*a6 + 0.44970998105242493*a7 - 0.44680890883509927*a8 + 0.3691002339464134*a9);

	double g0m2040DividedBym2 = (g0m2040DividedBym2LogZTerm * z2Logz + g0m2040DividedBym2Zm2Term) / meqHat2 + g0m2040DividedBym2ConstantTerm;
//	g0m2040DividedBym2=0.0;
//	g0m2040DividedBym2=g0m2040DividedBym2ConstantTerm;

	//---------------------------------------------------------------

	double g0m2020ConstantTerm = (-0.000011732877574172216 + 0.0008995982417255349*a + 0.7361757170685882*a2 + 8.855369493730423*a3 - 76.00475510029524*a4 + 261.7208147369771*a5 - 
     479.7229720410445*a6 + 487.6088140168412*a7 - 260.7802309619955*a8 + 57.58861537142237*a9)/
   (-0.0001445869624012341 + 0.03556951252765837*a + 3.2344009306364274*a2 + 12.514996006390815*a3 - 118.07971915906703*a4 + 347.32729764112844*a5 - 
     483.55235351080324*a6 + 310.7815908976544*a7 - 67.06626623263872*a8 - 5.19154598317511*a9);

	double g0m2020Z2Term = (-0.02907781069921067 - 13.745970067614193*a + 13997.84367889844*a10 - 325.2706635847172*a2 + 2604.8287335227815*a3 - 4049.1372422838062*a4 - 
     16667.365141623744*a5 + 83158.20010908133*a6 - 154605.3533428809*a7 + 147412.71289551136*a8 - 71512.66891762966*a9)/
   (0.9697656616566325 + 117.19185209632286*a + 7621.0452633237755*a10 + 842.3723825380051*a2 - 5926.237511293323*a3 + 14247.062719381864*a4 - 
     16353.128301113404*a5 + 21213.272784815457*a6 - 46076.07167476786*a7 + 59110.92666130391*a8 - 34796.81744513102*a9);

	double g0m2020Z2LogZTerm = (-0.0005534368375549849 - 0.3080130050685187*a + 630.3245252759478*a10 - 11.136413628061584*a2 + 34.63107256159966*a3 + 367.6376955018988*a4 - 
     2567.818354775122*a5 + 7178.293271719852*a6 - 10678.12443252271*a7 + 8758.211709132715*a8 - 3711.7083581405914*a9)/
   (0.009248889228618452 + 1.4181948328169438*a + 95.39051233960058*a10 + 21.503144137961534*a2 - 62.705284754794086*a3 - 83.17143447631213*a4 + 
     693.6443740085989*a5 - 1134.2586687505138*a6 + 513.1520937342799*a7 + 331.751797652403*a8 - 376.7067817672707*a9);

	double g0m2020 = g0m2020ConstantTerm + g0m2020Z2LogZTerm * z2Logz + g0m2020Z2Term * meqHat2;
//	g0m2020 = g0m2020ConstantTerm;

	//---------------------------------------------------------------

	double G_F1 = (0.5677352777359898 + 214.62760364687549*a + 6118.914112227032*a2 - 21609.531248796273*a3 + 20940.118061881778*a4 + 4120.31854449261*a5 - 18170.10263026395*a6 + 
     10045.925440196823*a7 - 1660.838410339361*a8)/
   (1.0013422651071762 + 216.42758743455585*a + 4808.207996489351*a2 - 13403.898599065215*a3 + 3036.0036716476625*a4 + 21220.315761322887*a5 - 
     24393.3243078871*a6 + 9629.478637585*a7 - 1114.218342270134*a8);

	double m2g1m2011 = G_F1 * meqHat2;

	//---------------------------------------------------------------

	double G_G1 = (0.000034723242182175396 + 0.012572522745974798*a + 0.427982420164622*a2 + 1.3191165093454609*a3 - 5.616153984227332*a4 + 4.28680523098434*a5 + 
     1.5082210795675917*a6 - 2.6464367629253513*a7 + 0.707946904809862*a8)/
   (0.00017101953273167538 + 0.030287072809487548*a + 0.5997230938576861*a2 + 0.3789662678043767*a3 - 5.429830325143101*a4 + 7.101821088278835*a5 - 
     1.968210179097073*a6 - 1.3532147373075232*a7 + 0.6402914643622666*a8);

	double m2g2m2000 = G_G1 * meqHat2;

	//---------------------------------------------------------------

	double G_H1 = (-6.155680471557033e-6 - 0.0437815902254004*a - 1.9202067785203083*a2 - 8.214545748547005*a3 + 16.851413984816478*a4 + 1.4197959008499166*a5 - 
     12.522032869988193*a6 + 4.396827516025351*a7)/
   (0.01623316424499292 + 0.9209657534005077*a + 6.978034999391149*a2 - 4.036784979260694*a3 - 20.42342324303473*a4 + 20.867517022584355*a5 - 
     2.9518500939454357*a6 - 1.3381486056401877*a7);

	double g1m2031 = G_H1;

	//---------------------------------------------------------------

	double G_L1 = (8.37681857229285e-6 - 0.037393624101204836*a - 0.08066673411547964*a2 + 10.588811403848602*a3 - 31.475314641340006*a4 + 33.88875878723781*a5 - 
     15.02360311716185*a6 + 2.139658098121864*a7)/
   (-0.203841931631807 + 0.12072040473285238*a + 48.38863379678098*a2 - 147.37851460276426*a3 + 153.83086096428207*a4 - 53.087554955269894*a5 - 
     6.884801658373992*a6 + 5.215153036261099*a7);

	double g2m2020 = G_L1;
	//======================================================================

	*beta_lPi = 3*g0m2040DividedBym2;
	*delta_lPi = 3/2*(g0m2020+g0m2040DividedBym2);
	*lambda_piPi = 3/4*(m2g0m2000+2*(g0m2020+1)+g0m2040DividedBym2);
	*beta_PiPi = 1-g0m2020ConstantTerm;
	*delta_PiPi = 0.5*(1+g0m2020ConstantTerm);
	*lambda_Pipi = 3*G_G1/5*(1/3-cs2);
}




