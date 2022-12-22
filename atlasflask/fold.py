import os


def trimmedData(seq,dms):
    dmsarr=dms.split(";")
    with open("./seq2",'w') as f1:
        f1.write(seq)

    with open("./dms2",'w') as f2:
        for i in range(len(dmsarr)):
            if (int(dmsarr[i])==0):
                dmsarr[i]=-999
            f2.write(str(i+1)+'\t'+str(dmsarr[i])+"\n")
    
    
seq='GTTGAGGGGAATGTTGTCTGGATCGAGGATATTATAGATATATACATGTGTATGTTAATGATTCAAGTGATCATAGAGAGTATCCTCGGACCAGGCTTCATCCCCCCCAAC'
dms='825;1791;2397;8097;1884;2360;2367;3658;1369;2006;913;4764;2114;5657;1677;1034;724;1323;1049;780;1314;1060;3442;1451;1221;1431;676;134;127;489;1018;1989;1260;7825;2154;1052;692;877;1483;1904;801;1318;2080;791;6082;2182;1558;1408;1580;1307;1461;1040;2258;2240;2181;3158;2585;845;1584;371;1924;567;970;1636;2378;2025;1807;736;638;687;78;548;420;719;1182;2805;1586;1950;371;419;663;16;3;3;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0'
trimmedData(seq,dms)
os.system('Fold.exe seq2 out2.ct -dms dms2')
os.system('ct2dot.exe out2.ct 1 out2.dot')