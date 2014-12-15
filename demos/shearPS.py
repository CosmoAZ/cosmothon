import sys
import math
import getopt
import cosmocalcs
import lensingClassFile
import numpy 
import matplotlib.pyplot as plt
from numpy import roll
from numpy import linspace
from numpy import logspace
from scipy.integrate import quad
from scipy.interpolate import interp1d
import powerspec

def main(argv):
    print ''
    print '=== Running nonlinearPS.py program ====\n'

    print argv

    #try:
    #opts, args = getopt.getopt(argv,"ho:")
    #except getopt.GetoptError as err: # if include option that's not there
    #    print "Program parameter not found\n"

    i = int(argv[3])
    j = int(argv[5])
    shearBins = int(argv[1])   

    f = open('Pk.txt', 'w')
    f.close()    
    fg = open('l.txt', 'w')   
    fg.close()

    """for opt, arg in opts:
        if opt == '-i':
            i = arg
        elif opt == '-j':
            j = arg
        elif opt in ("-n"):
            shearBins = arg
        elif opt in ("-h"):
            print "Run program with: python nonlinearPS.py -n <number of shear bins> -i <first bin to calculate cross-spectra spectrum> -j <second bin to calculate cross-spectra spectrum>"

    if i == '' or j == '' or shearBins == '':
	print "i, j, or shearBins not found\n"""

    print "Number of bins is: ", shearBins, " i is: ", i, " j is: ", j

    print "Initializing lc class\n"

    lc = lensingClassFile.lensingClass(0)

    print "Extracting redshifts from DLS data"

    F2redshifts = lc.getRedshiftsFromDLS()    	
    organizedListi, numberedListi = lc.organizeRedshifts(F2redshifts)
    sortedF2Redshifts = lc.sort_Redshift_Values(organizedListi)

    ns = 0.96  # spectral index
    DelRsq = 2.1e-9  # amplitude of primordial curvature fluctuations at the pivot scale

    #larray = linspace(0, 1000, 50)
    larray = logspace(math.log10(1), math.log10(1e5), 10000)

    pKij = []

    #Warraytt = [0.000295913, 0.000586919, 0.000873113, 0.001154586, 0.001431429, 0.001703731, 0.001971577, 0.002235053, 0.002494239, 0.002749218, 0.003000067, 0.003246863, 0.003489682, 0.003728597, 0.00396368 , 0.004195002, 0.00442263 , 0.004646634, 0.004867077, 0.005084025, 0.00529754 , 0.005507685, 0.005714519, 0.005918101, 0.00611849 , 0.006315742, 0.006509911, 0.006701053, 0.00688922 , 0.007074464, 0.007256836, 0.007436386, 0.007613163, 0.007787214, 0.007958585, 0.008127324, 0.008293475, 0.008457081, 0.008618186, 0.008776832, 0.00893306 , 0.009086911, 0.009238425, 0.00938764 , 0.009534595, 0.009679328, 0.009821874, 0.00996227 , 0.01010055 , 0.010236751, 0.010370905, 0.010503046, 0.010633206, 0.010761417, 0.010887712, 0.011012119, 0.011134671, 0.011255395, 0.011297328, 0.011336169, 0.011320878, 0.011302105, 0.011263476, 0.011221468, 0.01128943 , 0.011355324, 0.011281618, 0.011204851, 0.011239451, 0.01127189 , 0.011222967, 0.011171445, 0.011153462, 0.011133194, 0.011068391, 0.011001304, 0.010974857, 0.010946383, 0.010893488, 0.010838622, 0.010821445, 0.010802531, 0.010521606, 0.010241236, 0.010228996, 0.010215309, 0.010121096, 0.01002551 , 0.009998197, 0.009969566, 0.009902282, 0.00983375 , 0.009793138, 0.009751346, 0.009632994, 0.009513819, 0.009488328, 0.009461824, 0.009355087, 0.009247622, 0.00921389 , 0.009179259, 0.009107452, 0.009034869, 0.008973586, 0.008911551, 0.008788149, 0.008664474, 0.008486313, 0.008308749, 0.008182768, 0.008056744, 0.007995036, 0.007932822, 0.007835249, 0.007737464, 0.007644989, 0.007552314, 0.007474435, 0.007396288, 0.007323546, 0.00725054 , 0.007152876, 0.007055216, 0.006971444, 0.006887586, 0.006779605, 0.006671861, 0.006576421, 0.006481116, 0.006378491, 0.006276146, 0.006201256, 0.006126367, 0.006033466, 0.00594081 , 0.005859068, 0.005777477, 0.005692764, 0.005608277, 0.005544388, 0.005480526, 0.005384742, 0.005289421, 0.005243807, 0.005198119, 0.005144367, 0.005090621, 0.004984144, 0.004878443, 0.004819717, 0.004761095, 0.004700546, 0.004640142, 0.00455368 , 0.004467762, 0.004397494, 0.004327552, 0.00428386 , 0.004240206, 0.004190718, 0.004141339, 0.004087442, 0.004033719, 0.003950818, 0.003868571, 0.003775992, 0.003684328, 0.00362064 , 0.003557345, 0.00348318 , 0.003409631, 0.003315543, 0.003222593, 0.003191034, 0.003159528, 0.003083295, 0.003007851, 0.002965414, 0.002923176, 0.002878324, 0.002833719, 0.002774175, 0.002715156, 0.002666271, 0.002617737, 0.002558639, 0.002500121, 0.002447523, 0.002395399, 0.002337989, 0.002281191, 0.002235339, 0.002189882, 0.002145247, 0.002101005, 0.002052047, 0.0020036  , 0.001956929, 0.001910748, 0.001881382, 0.001852197, 0.001818893, 0.001785845, 0.001761792, 0.001737864, 0.001695539, 0.001653689, 0.001619578, 0.001585783, 0.001566919, 0.001548139, 0.001533879, 0.00151966 , 0.001489864, 0.001460332, 0.001438673, 0.00141715 , 0.001398229, 0.001379412, 0.001368522, 0.001357655, 0.001328641, 0.001299915, 0.001280768, 0.001261743, 0.001250906, 0.001240099, 0.001230419, 0.001220761, 0.001208863, 0.001197008, 0.001173141, 0.001149498, 0.001136786, 0.001124132, 0.001113348, 0.001102603, 0.001083036, 0.001063631, 0.001052727, 0.001041867, 0.001031317, 0.00102081 , 0.001006012, 0.000991312, 0.000977496, 0.000963767, 0.000942578, 0.000921614, 0.000892897, 0.000864625, 0.000855391, 0.0008462  , 0.000833638, 0.000821162, 0.000804399, 0.000787802, 0.000772356, 0.000757057, 0.000743332, 0.000729727, 0.000718761, 0.000707872, 0.000694623, 0.000681494, 0.000670852, 0.00066029 , 0.000652153, 0.000644062, 0.000629262, 0.000614631, 0.000604252, 0.000593958, 0.000587217, 0.000580511, 0.000567682, 0.000554994, 0.000542925, 0.000530987, 0.000516329, 0.000501874, 0.000494393, 0.000486965, 0.000480762, 0.000474597, 0.000442757, 0.000412025, 0.000406379, 0.000400771, 0.000395803, 0.000390864, 0.00038935 , 0.000387837, 0.000381691, 0.000375593, 0.000373359, 0.000371131, 0.000357618, 0.000344357, 0.000312869, 0.000282895, 0.000279936, 0.000276991, 0.000274962, 0.000272939, 0.000258394, 0.000244248, 0.000240369, 0.00023652 , 0.000235917, 0.000235314, 0.000217926, 0.000201208, 0.000197201, 0.000193234, 0.000184118, 0.000175224, 0.000174591, 0.00017396 , 0.000173201, 0.000172443, 0.000172183, 0.000171922, 0.000161934, 0.000152247, 0.00015177 , 0.000151293, 0.00015066 , 0.000150028, 0.000145723, 0.000141481, 0.000136192, 0.000131004, 0.000130346, 0.000129689, 0.000125469, 0.00012132 , 0.000118105, 0.000114934, 0.000113366, 0.000111809, 0.000107553, 0.000103382, 0.000102964, 0.000102548, 9.88062E-05, 9.51351E-05, 9.46362E-05, 9.41385E-05, 8.96393E-05, 8.52516E-05, 7.53495E-05, 6.60615E-05, 6.21952E-05, 5.84468E-05, 5.63407E-05, 5.42738E-05, 5.27104E-05, 5.11703E-05, 4.44793E-05, .82591E-05, 3.78342E-05, 3.74116E-05, 3.19856E-05, 2.69862E-05, 2.68388E-05, 2.66919E-05, 2.57949E-05, 2.49134E-05, 2.06624E-05, 1.68101E-05, 1.63942E-05, 1.59837E-05, 1.59004E-05, 1.58173E-05, 1.57043E-05, 1.55917E-05, 1.55533E-05, 1.5515E-05 , 1.54548E-05, 1.53947E-05, 1.48568E-05, 1.43287E-05, 1.4294E-05 , 1.42593E-05, 1.38513E-05, 1.34494E-05, 1.26638E-05, 1.19021E-05, 1.17943E-05, 1.16871E-05, 1.1153E-05 , 1.06315E-05, 9.76296E-06, 8.93174E-06, 8.88413E-06, 8.83665E-06, 7.82712E-06, 6.87919E-06, 6.56632E-06, 6.26086E-06, 6.2467E-06 , 6.23256E-06, 6.22597E-06, 6.21937E-06, 6.04572E-06, 5.8746E-06 , 5.86582E-06, 5.85703E-06, 5.84581E-06, 5.8346E-06 , 5.82704E-06, 5.81948E-06, 5.81555E-06, 5.81162E-06, 5.80406E-06, 5.79649E-06, 5.78651E-06, 5.77653E-06, 5.76776E-06, 5.759E-06  , 5.75264E-06, 5.74629E-06, 5.73393E-06, 5.72158E-06, 5.71283E-06, 5.70409E-06, 5.70014E-06, 5.69618E-06, 5.68505E-06, 5.67394E-06, 5.66878E-06, 5.66363E-06, 5.63589E-06, 5.60822E-06, 5.57708E-06, 5.54604E-06, 5.54095E-06, 5.53585E-06, 5.53076E-06, 5.52566E-06, 5.51938E-06, 5.51311E-06, 5.48695E-06, 5.46086E-06, 5.43949E-06, 5.41816E-06, 5.36336E-06, 5.30885E-06, 5.11623E-06, 4.92725E-06, 4.92144E-06, 4.91562E-06, 4.9043E-06 , 4.89298E-06, 4.85205E-06, 4.81131E-06, 4.58578E-06, 4.36578E-06, 4.29425E-06, 4.22335E-06, 4.21403E-06, 4.20472E-06, 4.17613E-06, 4.14764E-06, 4.13337E-06, 4.11913E-06, 4.11497E-06, 4.11081E-06, 4.10063E-06, 4.09047E-06, 4.04238E-06, 3.99459E-06, 3.97076E-06, 3.947E-06  , 3.94394E-06, 3.94088E-06, 3.93291E-06, 3.92495E-06, 3.89746E-06, 3.87007E-06, 3.86705E-06, 3.86402E-06, 3.861E-06  , 3.85797E-06, 3.84914E-06, 3.84031E-06, 3.83439E-06, 3.82848E-06, 3.82643E-06, 3.82438E-06, 3.81655E-06, 3.80873E-06, 3.80188E-06, 3.79504E-06, 3.79204E-06, 3.78904E-06, 3.78699E-06, 3.78495E-06, 3.7829E-06 , 3.78085E-06, 3.7788E-06 , 3.77676E-06, 3.77375E-06, 3.77075E-06, 3.76775E-06, 3.76474E-06, 3.76269E-06, 3.76064E-06, 3.75859E-06, 3.75654E-06, 3.75449E-06, 3.75243E-06, 3.75038E-06, 3.74833E-06, 3.74532E-06, 3.74232E-06, 3.74027E-06, 3.73821E-06, 3.73616E-06, 3.7341E-06 , 3.73205E-06, 3.72999E-06, 3.72794E-06, 3.72588E-06, 3.72382E-06, 3.72177E-06, 3.71971E-06, 3.71765E-06, 3.71559E-06, 3.71354E-06, 3.71148E-06, 3.70942E-06, 3.70736E-06, 3.7053E-06 , 3.70325E-06, 3.70119E-06, 3.69913E-06, 3.69707E-06, 3.69501E-06, 3.69296E-06, 3.6909E-06 , 3.68884E-06, 3.68678E-06, 3.68472E-06, 3.68267E-06, 3.68061E-06, 3.67855E-06, 3.6765E-06 , 3.67444E-06, 3.67238E-06, 3.67033E-06, 3.66827E-06, 3.66621E-06, 3.66416E-06, 3.6621E-06 , 3.66005E-06, 3.65799E-06, 3.65594E-06, 3.65388E-06, 3.65183E-06, 3.64978E-06, 3.64772E-06, 3.64567E-06, 3.64362E-06, 3.64157E-06, 3.63952E-06, 3.63746E-06, 3.63541E-06, 3.63336E-06, 3.63132E-06, 3.62927E-06, 3.62722E-06, 3.62425E-06, 3.62129E-06, 3.61924E-06, 3.6172E-06 , 3.61515E-06, 3.61311E-06, 3.61107E-06, 3.60902E-06, 3.60607E-06, 3.60312E-06, 3.60016E-06, 3.59721E-06, 3.59154E-06, 3.58587E-06, 3.58202E-06, 3.57817E-06, 3.57614E-06, 3.57412E-06, 3.57209E-06, 3.57006E-06, 3.56532E-06, 3.56059E-06, 3.55857E-06, 3.55655E-06, 3.55363E-06, 3.55071E-06, 3.5289E-06 , 3.50717E-06, 3.5016E-06 , 3.49603E-06, 3.49405E-06, 3.49206E-06, 3.47671E-06, 3.4614E-06 , 3.45766E-06, 3.45392E-06, 3.45018E-06, 3.44645E-06, 3.43388E-06, 3.42135E-06, 3.41764E-06, 3.41393E-06, 3.39267E-06, 3.37149E-06, 3.3687E-06 , 3.36591E-06, 3.35702E-06, 3.34814E-06, 3.33321E-06, 3.31832E-06, 3.31383E-06, 3.30935E-06, 3.30229E-06, 3.29524E-06, 3.28906E-06, 3.28289E-06, 3.269E-06  , 3.25516E-06, 3.24135E-06, 3.22757E-06, 3.21978E-06, 3.21201E-06, 3.1907E-06 , 3.16948E-06, 3.16599E-06, 3.1625E-06 , 3.15398E-06, 3.14547E-06, 3.11025E-06, 3.07525E-06, 3.06936E-06, 3.06348E-06, 3.03948E-06, 3.01558E-06, 2.99179E-06, 2.96811E-06, 2.95749E-06, 2.9469E-06 , 2.89771E-06, 2.84895E-06, 2.82829E-06, 2.80771E-06, 2.7637E-06 , 2.72006E-06, 2.70534E-06, 2.69067E-06, 2.6676E-06 , 2.64463E-06, 2.58539E-06, 2.52685E-06, 2.50902E-06, 2.49127E-06, 2.44121E-06, 2.39169E-06, 2.3744E-06 , 2.35718E-06, 2.34506E-06, 2.33298E-06, 2.2494E-06 , 2.16739E-06, 2.14896E-06, 2.13062E-06, 2.09066E-06, 2.0511E-06 , 2.0279E-06 , 2.00484E-06, 1.98851E-06, 1.97225E-06, 1.92354E-06, 1.87546E-06, 1.83687E-06, 1.7987E-06 , 1.76095E-06, 1.72362E-06, 1.70129E-06, 1.67911E-06, 1.65469E-06, 1.63046E-06, 1.61115E-06, 1.59197E-06, 1.5537E-06 , 1.51591E-06, 1.44654E-06, 1.37883E-06, 1.36443E-06, 1.35012E-06, 1.33159E-06, 1.31321E-06, 1.27233E-06, 1.23211E-06, 1.20885E-06, 1.18582E-06, 1.16904E-06, 1.15239E-06, 1.12256E-06, 1.09313E-06, 1.07322E-06, 1.05351E-06, 1.02693E-06, 1.0007E-06 , 9.82621E-07, 9.64713E-07, 9.45629E-07, 9.26746E-07, 9.01027E-07, 8.75684E-07, 8.51996E-07, 8.28645E-07, 8.11031E-07, 7.93617E-07, 7.69905E-07, 7.46566E-07, 7.35027E-07, 7.23584E-07, 6.9483E-07 , 6.66676E-07, 6.50964E-07, 6.35449E-07, 6.00695E-07, 5.66938E-07, 5.5902E-07 , 5.51163E-07, 5.39979E-07, 5.28916E-07, 5.19298E-07, 5.09773E-07, 4.81664E-07, 4.54367E-07, 4.3482E-07 , 4.15713E-07, 4.06356E-07, 3.97112E-07, 3.84272E-07, 3.7165E-07 , 3.69209E-07, 3.66777E-07, 3.49027E-07, 3.31727E-07, 3.27073E-07, 3.22454E-07, 3.20456E-07, 3.18466E-07, 3.16483E-07, 3.14508E-07, 2.9689E-07 , 2.7979E-07 , 2.76741E-07, 2.7371E-07 , 2.69511E-07, 2.65346E-07, 2.63086E-07, 2.60837E-07, 2.59528E-07, 2.58223E-07, 2.45043E-07, 2.32215E-07, 2.25328E-07, 2.1855E-07 , 2.18002E-07, 2.17456E-07, 2.15006E-07, 2.12572E-07, 2.06209E-07, 1.99945E-07, 1.89997E-07, 1.80308E-07, 1.79244E-07, 1.78184E-07, 1.74841E-07, 1.71532E-07, 1.69936E-07, 1.68349E-07, 1.64553E-07, 1.60803E-07, 1.44941E-07, 1.2991E-07 , 1.29186E-07, 1.28464E-07, 1.26613E-07, 1.24777E-07, 1.23115E-07, 1.21464E-07, 1.20296E-07, 1.19134E-07, 1.08991E-07, 9.93049E-08, 9.79726E-08, 9.665E-08  , 9.61759E-08, 9.57032E-08, 9.52319E-08, 9.4762E-08 , 8.72259E-08, 8.00059E-08, 7.88182E-08, 7.76399E-08, 7.49828E-08, 7.23733E-08, 6.80374E-08, 6.38378E-08, 6.34619E-08, 6.30874E-08, 5.97083E-08, 5.64241E-08, 5.48024E-08, 5.32052E-08, 4.94099E-08, 4.57569E-08, 4.52527E-08, 4.47515E-08, 3.79759E-08, 3.17594E-08, 3.0554E-08 , 2.93726E-08, 2.79895E-08, 2.66405E-08, 2.63342E-08, 2.60297E-08, 2.57992E-08, 2.55697E-08, 2.30437E-08, 2.06503E-08, 2.01914E-08, 1.97379E-08, 1.97271E-08, 1.97164E-08, 1.95178E-08, 1.93203E-08, 1.81481E-08, 1.70132E-08, 1.70039E-08, 1.69947E-08, .66382E-08, 1.62857E-08, 1.48308E-08, 1.34447E-08, 1.32827E-08, 1.31219E-08, 1.30129E-08, 1.29044E-08, 1.2496E-08 , 1.20943E-08, 1.14599E-08, 1.0843E-08 , 1.0561E-08 , 1.02829E-08, 1.02773E-08, 1.02718E-08, 1.00423E-08, 9.81557E-09, 9.63509E-09, 9.45637E-09, 9.40815E-09, 9.36008E-09, 9.2694E-09, 9.17921E-09, 9.04728E-09, 8.91637E-09, 8.33546E-09, 7.77441E-09, 7.46099E-09, 7.15418E-09, 7.15033E-09, 7.14648E-09, 7.10529E-09, 7.06424E-09, 6.87586E-09, 6.69013E-09, 6.50704E-09, 6.32659E-09, 5.90868E-09, 5.50525E-09, 5.27543E-09, 5.05062E-09, 4.64855E-09, 4.26335E-09, 4.14669E-09, 4.0317E-09, 3.97377E-09, 3.91629E-09, 3.72362E-09, 3.53591E-09, 3.48186E-09, 3.42825E-09, 3.34957E-09, 3.27184E-09, 3.19507E-09, 3.11925E-09, 2.76024E-09, 2.42335E-09, 2.29404E-09, 2.16835E-09, 2.1672E-09, 2.16604E-09, 1.96547E-09, 1.77475E-09, 1.70066E-09, 1.6282E-09, 1.62734E-09, 1.62647E-09, 1.43705E-09, 1.25945E-09, 1.16727E-09, 1.07863E-09, 9.9354E-10, 9.11987E-10, 8.59429E-10, 8.08458E-10, 7.83359E-10, 7.58669E-10, 7.2259E-10, 6.87409E-10, 5.98493E-10, 5.15776E-10, 4.86192E-10, .57497E-10, 4.47972E-10, 4.38552E-10, 4.20248E-10, 4.02343E-10, 3.93437E-10, 3.84635E-10, 3.42904E-10, 3.0359E-10, 2.95892E-10, 2.88297E-10, 2.88145E-10, 2.87994E-10, 2.8051E-10, 2.73128E-10, 2.45006E-10, 2.18424E-10, 2.11936E-10, 2.05549E-10, 1.99264E-10, 1.93078E-10, 1.7531E-10, 1.58408E-10, 1.47594E-10, 1.37168E-10, 1.22287E-10, 1.08267E-10, 1.0821E-10, 1.08154E-10, 1.0364E-10, 9.92253E-11, 8.26759E-11, 6.76428E-11, 6.76076E-11, 6.75724E-11, 6.06152E-11, 5.40392E-11, 4.78437E-11, 4.20283E-11, 3.65923E-11, 3.15352E-11, 2.68563E-11, 2.25551E-11, 2.05405E-11, 1.86212E-11, 1.86116E-11, 1.86019E-11, 1.3433E-11, 9.10553E-12, 9.10082E-12, 9.09612E-12, 6.67941E-12, 4.63609E-12, 3.75329E-12, 2.96404E-12, 2.96251E-12, 2.96098E-12, 2.26583E-12, 1.66384E-12, 1.15485E-12, 7.38721E-13, 4.15317E-13, 1.8449E-13, 4.60989E-14, 0, 0, 0, 0, 0]
    Warraytt =[0.001479565,0.002934595,0.004365565,0.00577293,0.007157145,0.008518655,0.009857885,0.011175265,0.012471195,0.01374609,0.015000335,0.016234315,0.01744841,0.018642985,0.0198184,0.02097501,0.02211315,0.02323317,0.024335385,0.025420125,0.0264877,0.027538425,0.028572595,0.029590505,0.03059245,0.03157871,0.032549555,0.033505265,0.0344461,0.03537232,0.03628418,0.03718193,0.038065815,0.03893607,0.039792925,0.04063662,0.041467375,0.042285405,0.04309093,0.04388416,0.0446653,0.045434555,0.046192125,0.0469382,0.047672975,0.04839664,0.04910937,0.04981135,0.05050275,0.051183755,0.051854525,0.05251523,0.05316603,0.053807085,0.05443856,0.055060595,0.055673355,0.056276975,0.05648664,0.056680845,0.05660439,0.056510525,0.05631738,0.05610734,0.05644715,0.05677662,0.05640809,0.056024255,0.056197255,0.05635945,0.056114835,0.055857225,0.05576731,0.05566597,0.055341955,0.05500652,0.054874285,0.054731915,0.05446744,0.05419311,0.054107225,0.054012655,0.05260803,0.05120618,0.05114498,0.051076545,0.05060548,0.05012755,0.049990985,0.04984783,0.04951141,0.04916875,0.04896569,0.04875673,0.04816497,0.047569095,0.04744164,0.04730912,0.046775435,0.04623811,0.04606945,0.045896295,0.04553726,0.045174345,0.04486793,0.044557755,0.043940745,0.04332237,0.042431565,0.041543745,0.04091384,0.04028372,0.03997518,0.03966411,0.039176245,0.03868732,0.038224945,0.03776157,0.037372175,0.03698144,0.03661773,0.0362527,0.03576438,0.03527608,0.03485722,0.03443793,0.033898025,0.033359305,0.032882105,0.03240558,0.031892455,0.03138073,0.03100628,0.030631835,0.03016733,0.02970405,0.02929534,0.028887385,0.02846382,0.028041385,0.02772194,0.02740263,0.02692371,0.026447105,0.026219035,0.025990595,0.025721835,0.025453105,0.02492072,0.024392215,0.024098585,0.023805475,0.02350273,0.02320071,0.0227684,0.02233881,0.02198747,0.02163776,0.0214193,0.02120103,0.02095359,0.020706695,0.02043721,0.020168595,0.01975409,0.019342855,0.01887996,0.01842164,0.0181032,0.017786725,0.0174159,0.017048155,0.016577715,0.016112965,0.01595517,0.01579764,0.015416475,0.015039255,0.01482707,0.01461588,0.01439162,0.014168595,0.013870875,0.01357578,0.013331355,0.013088685,0.012793195,0.012500605,0.012237615,0.011976995,0.011689945,0.011405955,0.011176695,0.01094941,0.010726235,0.010505025,0.010260235,0.010018,0.009784645,0.00955374,0.00940691,0.009260985,0.009094465,0.008929225,0.00880896,0.00868932,0.008477695,0.008268445,0.00809789,0.007928915,0.007834595,0.007740695,0.007669395,0.0075983,0.00744932,0.00730166,0.007193365,0.00708575,0.006991145,0.00689706,0.00684261,0.006788275,0.006643205,0.006499575,0.00640384,0.006308715,0.00625453,0.006200495,0.006152095,0.006103805,0.006044315,0.00598504,0.005865705,0.00574749,0.00568393,0.00562066,0.00556674,0.005513015,0.00541518,0.005318155,0.005263635,0.005209335,0.005156585,0.00510405,0.00503006,0.00495656,0.00488748,0.004818835,0.00471289,0.00460807,0.004464485,0.004323125,0.004276955,0.004231,0.00416819,0.00410581,0.004021995,0.00393901,0.00386178,0.003785285,0.00371666,0.003648635,0.003593805,0.00353936,0.003473115,0.00340747,0.00335426,0.00330145,0.003260765,0.00322031,0.00314631,0.003073155,0.00302126,0.00296979,0.002936085,0.002902555,0.00283841,0.00277497,0.002714625,0.002654935,0.002581645,0.00250937,0.002471965,0.002434825,0.00240381,0.002372985,0.002213785,0.002060125,0.002031895,0.002003855,0.001979015,0.00195432,0.00194675,0.001939185,0.001908455,0.001877965,0.001866795,0.001855655,0.00178809,0.001721785,0.001564345,0.001414475,0.00139968,0.001384955,0.00137481,0.001364695,0.00129197,0.00122124,0.001201845,0.0011826,0.001179585,0.00117657,0.00108963,0.00100604,0.000986005,0.00096617,0.00092059,0.00087612,0.000872955,0.0008698,0.000866005,0.000862215,0.000860915,0.00085961,0.00080967,0.000761235,0.00075885,0.000756465,0.0007533,0.00075014,0.000728615,0.000707405,0.00068096,0.00065502,0.00065173,0.000648445,0.000627345,0.0006066,0.000590525,0.00057467,0.00056683,0.000559045,0.000537765,0.00051691,0.00051482,0.00051274,0.000494031,0.000475676,0.000473181,0.000470693,0.000448197,0.000426258,0.000376748,0.000330308,0.000310976,0.000292234,0.000281704,0.000271369,0.000263552,0.000255852,0.000222397,4.12955E-05,0.000189171,0.000187058,0.000159928,0.000134931,0.000134194,0.00013346,0.000128975,0.000124567,0.000103312,8.40505E-05,0.000081971,7.99185E-05,0.000079502,7.90865E-05,7.85215E-05,7.79585E-05,7.77665E-05,0.000077575,0.000077274,7.69735E-05,0.000074284,7.16435E-05,0.00007147,7.12965E-05,6.92565E-05,0.000067247,0.000063319,5.95105E-05,5.89715E-05,5.84355E-05,0.000055765,5.31575E-05,4.88148E-05,4.46587E-05,4.44207E-05,4.41833E-05,3.91356E-05,3.4396E-05,3.28316E-05,3.13043E-05,3.12335E-05,3.11628E-05,3.11299E-05,3.10969E-05,3.02286E-05,0.000029373,2.93291E-05,2.92852E-05,2.92291E-05,0.000029173,2.91352E-05,2.90974E-05,2.90778E-05,2.90581E-05,2.90203E-05,2.89825E-05,2.89326E-05,2.88827E-05,2.88388E-05,0.000028795,2.87632E-05,2.87315E-05,2.86697E-05,2.86079E-05,2.85642E-05,2.85205E-05,2.85007E-05,2.84809E-05,2.84253E-05,2.83697E-05,2.83439E-05,2.83182E-05,2.81795E-05,2.80411E-05,2.78854E-05,2.77302E-05,2.77048E-05,2.76793E-05,2.76538E-05,2.76283E-05,2.75969E-05,2.75656E-05,2.74348E-05,2.73043E-05,2.71975E-05,2.70908E-05,2.68168E-05,2.65443E-05,2.55812E-05,2.46363E-05,2.46072E-05,2.45781E-05,2.45215E-05,2.44649E-05,2.42603E-05,2.40566E-05,2.29289E-05,2.18289E-05,2.14713E-05,2.11168E-05,2.10702E-05,2.10236E-05,2.08807E-05,2.07382E-05,2.06669E-05,2.05957E-05,2.05749E-05,2.05541E-05,2.05032E-05,2.04524E-05,2.02119E-05,1.9973E-05,1.98538E-05,0.000019735,1.97197E-05,1.97044E-05,1.96646E-05,1.96248E-05,1.94873E-05,1.93504E-05,1.93353E-05,1.93201E-05,0.000019305,1.92899E-05,1.92457E-05,1.92016E-05,1.9172E-05,1.91424E-05,1.91322E-05,1.91219E-05,1.90828E-05,1.90437E-05,1.90094E-05,1.89752E-05,1.89602E-05,1.89452E-05,1.8935E-05,1.89248E-05,1.89145E-05,1.89043E-05,0.000018894,1.88838E-05,1.88688E-05,1.88538E-05,1.88388E-05,1.88237E-05,1.88135E-05,1.88032E-05,1.8793E-05,1.87827E-05,1.87725E-05,1.87622E-05,1.87519E-05,1.87417E-05,1.87266E-05,1.87116E-05,1.87014E-05,1.86911E-05,1.86808E-05,1.86705E-05,1.86603E-05,1.865E-05,1.86397E-05,1.86294E-05,1.86191E-05,1.86089E-05,1.85986E-05,1.85883E-05,1.8578E-05,1.85677E-05,1.85574E-05,1.85471E-05,1.85368E-05,1.85265E-05,1.85163E-05,1.8506E-05,1.84957E-05,1.84854E-05,1.84751E-05,1.84648E-05,1.84545E-05,1.84442E-05,1.84339E-05,1.84236E-05,1.84134E-05,1.84031E-05,1.83928E-05,1.83825E-05,1.83722E-05,1.83619E-05,1.83517E-05,1.83414E-05,1.83311E-05,1.83208E-05,1.83105E-05,1.83003E-05,1.829E-05,1.82797E-05,1.82694E-05,1.82592E-05,1.82489E-05,1.82386E-05,1.82284E-05,1.82181E-05,1.82079E-05,1.81976E-05,1.81873E-05,1.81771E-05,1.81668E-05,1.81566E-05,1.81464E-05,1.81361E-05,1.81213E-05,1.81065E-05,1.80962E-05,0.000018086,1.80758E-05,1.80656E-05,1.80554E-05,1.80451E-05,1.80304E-05,1.80156E-05,1.80008E-05,1.79861E-05,1.79577E-05,1.79294E-05,1.79101E-05,1.78909E-05,1.78807E-05,1.78706E-05,1.78605E-05,1.78503E-05,1.78266E-05,1.7803E-05,1.77929E-05,1.77828E-05,1.77682E-05,1.77536E-05,1.76445E-05,1.75359E-05,0.000017508,1.74802E-05,1.74703E-05,1.74603E-05,1.73836E-05,0.000017307,1.72883E-05,1.72696E-05,1.72509E-05,1.72323E-05,1.71694E-05,1.71068E-05,1.70882E-05,1.70697E-05,1.69634E-05,1.68575E-05,1.68435E-05,1.68296E-05,1.67851E-05,1.67407E-05,1.66661E-05,1.65916E-05,1.65692E-05,1.65468E-05,1.65115E-05,1.64762E-05,1.64453E-05,1.64145E-05,0.000016345,1.62758E-05,1.62068E-05,1.61379E-05,1.60989E-05,1.60601E-05,1.59535E-05,1.58474E-05,1.583E-05,1.58125E-05,1.57699E-05,1.57274E-05,1.55513E-05,1.53763E-05,1.53468E-05,1.53174E-05,1.51974E-05,1.50779E-05,1.4959E-05,1.48406E-05,1.47875E-05,1.47345E-05,1.44886E-05,1.42448E-05,1.41415E-05,1.40386E-05,1.38185E-05,1.36003E-05,1.35267E-05,1.34534E-05,0.000013338,1.32232E-05,1.2927E-05,1.26343E-05,1.25451E-05,1.24564E-05,1.22061E-05,1.19585E-05,0.000011872,1.17859E-05,1.17253E-05,1.16649E-05,0.000011247,1.0837E-05,1.07448E-05,1.06531E-05,1.04533E-05,1.02555E-05,1.01395E-05,1.00242E-05,9.94255E-06,9.86125E-06,9.6177E-06,9.3773E-06,9.18435E-06,8.9935E-06,8.80475E-06,8.6181E-06,8.50645E-06,8.39555E-06,8.27345E-06,8.1523E-06,8.05575E-06,7.95985E-06,7.7685E-06,7.57955E-06,7.2327E-06,6.89415E-06,6.82215E-06,6.7506E-06,6.65795E-06,6.56605E-06,6.36165E-06,6.16055E-06,6.04425E-06,5.9291E-06,5.8452E-06,5.76195E-06,5.6128E-06,5.46565E-06,5.3661E-06,5.26755E-06,5.13465E-06,5.0035E-06,4.91311E-06,4.82357E-06,4.72815E-06,4.63373E-06,4.50514E-06,4.37842E-06,4.25998E-06,4.14323E-06,4.05516E-06,3.96809E-06,3.84953E-06,3.73283E-06,3.67514E-06,3.61792E-06,3.47415E-06,3.33338E-06,3.25482E-06,3.17725E-06,3.00348E-06,2.83469E-06,2.7951E-06,2.75582E-06,2.6999E-06,2.64458E-06,2.59649E-06,2.54887E-06,2.40832E-06,2.27184E-06,2.1741E-06,2.07857E-06,2.03178E-06,1.98556E-06,1.92136E-06,1.85825E-06,1.84605E-06,1.83389E-06,1.74514E-06,1.65864E-06,1.63537E-06,1.61227E-06,1.60228E-06,1.59233E-06,1.58242E-06,1.57254E-06,1.48445E-06,1.39895E-06,1.38371E-06,1.36855E-06,1.34756E-06,1.32673E-06,1.31543E-06,1.30419E-06,1.29764E-06,1.29112E-06,1.22522E-06,1.16108E-06,1.12664E-06,1.09275E-06,1.09001E-06,1.08728E-06,1.07503E-06,1.06286E-06,1.03105E-06,9.99725E-07,9.49985E-07,9.0154E-07,8.9622E-07,8.9092E-07,8.74205E-07,8.5766E-07,8.4968E-07,8.41745E-07,8.22765E-07,8.04015E-07,7.24705E-07,6.4955E-07,6.4593E-07,6.4232E-07,6.33065E-07,6.23885E-07,6.15575E-07,6.0732E-07,6.0148E-07,5.9567E-07,5.44955E-07,4.96525E-07,4.89863E-07,4.8325E-07,4.8088E-07,4.78516E-07,4.7616E-07,4.7381E-07,4.3613E-07,4.0003E-07,3.94091E-07,3.882E-07,3.74914E-07,3.61867E-07,3.40187E-07,3.19189E-07,3.1731E-07,3.15437E-07,2.98542E-07,2.82121E-07,2.74012E-07,2.66026E-07,2.4705E-07,2.28785E-07,2.26264E-07,2.23758E-07,1.8988E-07,1.58797E-07,1.5277E-07,1.46863E-07,1.39948E-07,1.33203E-07,1.31671E-07,1.30149E-07,1.28996E-07,1.27849E-07,1.15219E-07,1.03252E-07,1.00957E-07,9.86895E-08,9.86355E-08,9.8582E-08,9.7589E-08,9.66015E-08,9.07405E-08,8.5066E-08,8.50195E-08,8.49735E-08,3.3191E-08,8.14285E-08,7.4154E-08,6.72235E-08,6.64135E-08,6.56095E-08,6.50645E-08,6.4522E-08,6.248E-08,6.04715E-08,5.72995E-08,5.4215E-08,5.2805E-08,5.14145E-08,5.13865E-08,5.1359E-08,5.02115E-08,4.90779E-08,4.81755E-08,4.72819E-08,4.70408E-08,4.68004E-08,4.6347E-08,4.58961E-08,4.52364E-08,4.45819E-08,4.16773E-08,3.88721E-08,3.7305E-08,3.57709E-08,3.57517E-08,3.57324E-08,3.55265E-08,3.53212E-08,3.43793E-08,3.34507E-08,3.25352E-08,3.1633E-08,2.95434E-08,2.75263E-08,2.63772E-08,2.52531E-08,2.32428E-08,2.13168E-08,2.07335E-08,2.01585E-08,1.98689E-08,1.95815E-08,1.86181E-08,1.76796E-08,1.74093E-08,1.71413E-08,1.67479E-08,1.63592E-08,1.59754E-08,1.55963E-08,1.38012E-08,1.21168E-08,1.14702E-08,1.08418E-08,1.0836E-08,1.08302E-08,9.82735E-09,8.87375E-09,8.5033E-09,8.141E-09,8.1367E-09,8.13235E-09,7.18525E-09,6.29725E-09,5.83635E-09,5.39315E-09,4.9677E-09,4.55994E-09,4.29715E-09,4.04229E-09,3.9168E-09,3.79335E-09,3.61295E-09,3.43705E-09,2.99247E-09,2.57888E-09,2.43096E-09,2.87485E-10,2.23986E-09,2.19276E-09,2.10124E-09,2.01172E-09,1.96719E-09,1.92318E-09,1.71452E-09,1.51795E-09,1.47946E-09,1.44149E-09,1.44073E-09,1.43997E-09,1.40255E-09,1.36564E-09,1.22503E-09,1.09212E-09,1.05968E-09,1.02775E-09,9.9632E-10,9.6539E-10,8.7655E-10,7.9204E-10,7.3797E-10,6.8584E-10,6.11435E-10,5.41335E-10,5.4105E-10,5.4077E-10,5.182E-10,4.96127E-10,4.1338E-10,3.38214E-10,3.38038E-10,3.37862E-10,3.03076E-10,2.70196E-10,2.39219E-10,2.10142E-10,1.82962E-10,1.57676E-10,1.34282E-10,1.12776E-10,1.02703E-10,9.3106E-11,9.3058E-11,9.30095E-11,6.7165E-11,4.55277E-11,4.55041E-11,4.54806E-11,3.33971E-11,2.31805E-11,1.87665E-11,1.48202E-11,1.48126E-11,1.48049E-11,1.13292E-11,8.3192E-12,5.77425E-12,3.69361E-12,2.07659E-12,9.2245E-13,2.30495E-13,0,0,0,0,0]

    #Warraytt = Warraytt[0:3:len(Warraytt)]

    f_x = interp1d(linspace(0, max(sortedF2Redshifts), len(Warraytt)), Warraytt)

    redshifts = linspace(0, max(sortedF2Redshifts), len(sortedF2Redshifts))
    Warraytt = f_x(redshifts)

    f = open('stuffyear.txt', 'w')
    f.write(str(larray))
    f.close()

    #Warraytt = []
    Warrayi = []
    Warrayj = []

    if i > shearBins - 1 or j > shearBins - 1 or i < 0 or j < 0:

	print "i or j (bin number) is greater than the number of Shear bins: ", shearBins
	return

    shearbinarray = linspace(0, max(sortedF2Redshifts), shearBins+1)

    if shearBins == 1:
	zttiarray = F2redshifts

    	#organizedListi, numberedListi = lc.organizeRedshifts(zttiarray)
	#sortedF2Redshiftsi = lc.sort_Redshift_Values(organizedListi)

    else:
	zttiarray = []
	zttjarray = []

	for ztt in F2redshifts:
		if ztt >= shearbinarray[i] and ztt < shearbinarray[i+1]:
			zttiarray.append(ztt)

		if ztt >= shearbinarray[j] and ztt < shearbinarray[j+1]:
			zttjarray.append(ztt)

	organizedListi, numberedListi = lc.organizeRedshifts(zttiarray)
	sortedF2Redshiftsi = lc.sort_Redshift_Values(organizedListi)
	
	organizedListj, numberedListj = lc.organizeRedshifts(zttjarray)
	sortedF2Redshiftsj = lc.sort_Redshift_Values(organizedListj)

    k = 0

    """if shearBins == 1:
	for z in sortedF2Redshifts[1:len(sortedF2Redshifts) - 2]:
		g = open('lists.txt', 'a')
		g.write('orgListi: ' + str(organizedListi[k]) + '\t' + 'numListi: ' + str(numberedListi[k]) + '\n')
    		g.close()
		k += 1
    else:
	for z in sortedF2Redshifts[1:len(sortedF2Redshifts) - 1]:
		g = open('lists.txt', 'a')
		g.write('orgListi: ' + str(organizedListi[k]) + '\t' + 'orgListj: ' + str(organizedListj[k]) + '\t' + 'numListi: ' + str(numberedListi[k]) + '\t' +  'numListj: ' + str(numberedListj[k]) + '\n')
		k += 1
	    	g.close()"""

    Wvali = 0
    tempWvali = 0
    Wvalj = 0
    tempWvalj = 0

    """
    for z in sortedF2Redshifts[1:len(sortedF2Redshifts) - 1]:

    	print "Determining lensing weight function at z: ", z

	if shearBins == 1:

    		lc = lensingClassFile.lensingClass(0)
		Wvali, tempWvali = lc.lensingWeightFunction(z, sortedF2Redshiftsi, organizedListi, numberedListi)
		Warraytt.append(Wvali)	

	else:
	
		lc = lensingClassFile.lensingClass(0)
		Wvali, tempWvali = lc.lensingWeightFunction(z, sortedF2Redshiftsi, organizedListi, numberedListi)
		lc = lensingClassFile.lensingClass(0)
		Wvalj, tempWvalj = lc.lensingWeightFunction(z, sortedF2Redshiftsj, organizedListj, numberedListj)

		Warrayi.append(Wvali)
		Warrayj.append(Wvalj)
    """

    k = 0

    for z in sortedF2Redshifts[1:len(sortedF2Redshifts) - 1]:
        f = open('shearW.txt', 'a')
	if shearBins == 1:
		f.write('n: ' + str(numberedListi) + '\t' + 'W: ' + str(Warraytt[k]) + '\n')
	else:
		f.write('n_i: ' + str(numberedListi[k]) + '\t' + 'n_j: ' + str(numberedListj[k]) + '\t' + 'W_i: ' + str(Warrayi[k]) + '\t' + 'W_j: ' + str(Warrayj[k]) + '\n')
	f.close()
	k += 1

    cosmoCalcs = cosmocalcs.cosmologyCalculator(lc.h, lc.omegams[0], lc.omegams[1], lc.wX, lc.wXa)

    tf = powerspec.transferFunction(lc.omegams[0]-0.05, lc.omegams[2], lc.h)

    gro = powerspec.growthFunctionLCDM(lc.omegams[0], lc.omegams[1])

    ps = powerspec.powerSpectrum(tf, gro, ns, DelRsq)

    for l in larray:

	pKij_l = 0
	ii = 0

    	for z in sortedF2Redshifts[1:len(sortedF2Redshifts) - 1]:

		#print "Determining cross-spectra weight function at z: ", z

		#cosmoCalcs = cosmocalcs.cosmologyCalculator(lc.h, lc.omegams[0], lc.omegams[1], lc.wX, lc.wXa)

		cosmoCalcs.setEmissionRedShift(z)

		#D_A = (1 + z) * cosmoCalcs.AngularDiameterDistance()
		D_AMpc = (1 + z) * cosmoCalcs.AngularDiameterDistanceMpc()

		Pk = ps.getPk(l/D_AMpc, z) # in (Mpc/h)^3 units

		if shearBins == 1:
			pKij_l += Warraytt[ii]*Warraytt[ii]*Pk/(D_AMpc*D_AMpc*100*lc.h*pow(lc.omegams[0]*pow(1 + z, 3) + (1 - lc.omegams[1] - lc.omegams[2])*pow(1 + z, 2) + lc.omegams[2], 0.5))*(sortedF2Redshifts[1] - sortedF2Redshifts[0])
		else:
			pKij_l += Warrayi[ii]*Warrayj[ii]*Pk/(D_AMpc*D_AMpc*100*lc.h*pow(lc.omegams[0]*pow(1 + z, 3) + (1 - lc.omegams[1] - lc.omegams[2])*pow(1 + z, 2) + lc.omegams[2], 0.5))*(sortedF2Redshifts[1] - sortedF2Redshifts[0])

		ii += 1

	pKij.append(pKij_l)
	
        f = open('Pk.txt', 'a')
	fg = open('l.txt', 'a')
	#f.write('pKif_l: ' + str(pKij_l) + '\t' + 'l: ' + str(l) + '\t' + 'l(l+1): ' + str(l*(l + 1)) + '\n')
	f.write(str(pKij_l) + '\n')
	fg.write(str(l) + '\n')
	f.close()
	fg.close()


    """pKij_l_1 = []
    pKij_l_1_1bins = []
    i = 0

    for lt in larray:
	pKij_l_1.append(lt*(lt + 1)*pKij[i])
	pKij_l_1_1bins.append(lt*(lt + 1)*pKij[i]/(2.0*math.pi))
	i += 1

    plt.plot(larray, pKij, 'r-')
    #plt.xscale('symlog')
    #plt.yscale('symlog')
    plt.semilogy()
    plt.semilogx()
    plt.xlabel('l')
    plt.ylabel('$P_K^{' + str(i) + str(j) + '}$')
    plt.title('Shear Power Spectrum (' + str(shearBins) + ' bins) vs. Multipole Number')
    plt.axis([min(larray), max(larray), min(pKij), max(pKij)])
    plt.grid(True)
    plt.savefig('pKij_l.png')
    plt.close()

    plt.plot(larray, pKij_l_1, 'r-')
    #plt.xscale('symlog')
    #plt.yscale('symlog')
    plt.semilogy()
    plt.semilogx()
    plt.xlabel('l (multipole number)')
    plt.ylabel('$l(l + 1)/2\pi P_K^{' + str(i) + str(j) + '}$')
    plt.title('Shear Power Spectrum (' + str(shearBins) + ' bins) vs. Multipole Number')
    plt.axis([min(larray), max(larray), min(pKij_l_1), max(pKij_l_1)])
    plt.grid(True)
    plt.savefig('pKij_l_1.png')
    plt.close()"""


if __name__ == "__main__":
    main(sys.argv[1:])
    
