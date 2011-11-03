using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Runtime.InteropServices;
using System.Xml;
using System.Text.RegularExpressions;
using System.Net;
using System.Reflection;
using System.ComponentModel;
using pwiz.CLI.analysis;
using pwiz.CLI.data;
using pwiz.CLI.msdata;
using pwiz.CLI.proteome;
using SourceList = System.Collections.Generic.Set<IDPicker.SourceInfo>;



namespace IDPReader
{
    class IdpReader
    {
        //get pepSequence,source,scanID, write into a csv file. 
        //have tons of information: theretical ion intensity, pep info, chargeLabel...
        public static void idpReader(string idpXMLFile, string mzMLFile, double TicCutoffPercentage, int z, int model)
        {
            //get the path and filename of output csv file:
            string fileName = Path.GetFileNameWithoutExtension(idpXMLFile);
            string filePath = Path.GetDirectoryName(idpXMLFile);
            string csvFile = Path.Combine(filePath, fileName) + "_" + z.ToString() + ".csv";
            

            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, idpXMLFile);

            using (StreamWriter file = new StreamWriter(csvFile))
            {
                //TODO idOrIndex + "," + pepSequence + "," + pepBond + "," + AA + "," + bIons + "," + yIons;
                file.WriteLine("nativeID,pepSequence,bond,b1,b2,y1,y2");
                
                MSDataFile foo = new MSDataFile(mzMLFile);
                SpectrumList sl = foo.run.spectrumList;
                                    
                foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
                    foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
                        foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
                        {

                            IDPicker.ResultInstance ri = sItr.Value.results[1];
                            IDPicker.VariantInfo vi = ri.info.peptides.Min;
                            
                            string ss = vi.ToString() + "," + sItr.Value.id.source.name + "," + sItr.Value.nativeID;
                            bool boolCharge = sItr.Value.id.charge.Equals(z);
                            if (boolCharge)
                            {
                                string rawPepSequence = vi.ToString();
                                string pepSequence = vi.peptide.sequence;
                                int len = pepSequence.Length;
                                
                                // Look up the index with nativeID
                                object idOrIndex = null;
                                if( sItr.Value.nativeID != null && sItr.Value.nativeID.Length > 0 )
                                    idOrIndex = sItr.Value.nativeID;
                                int spectrumIndex = sl.find(idOrIndex as string);
                                // Trust the local index, if the nativeID lookup fails
                                if( spectrumIndex >= sl.size() )
                                    spectrumIndex = sItr.Value.id.index;
                                // Bail of the loca index is larger than the spectrum list size
                                if( spectrumIndex >= sl.size() )
                                    throw new Exception( "Can't find spectrum associated with the index." );

                                //get base peak and TIC and converted to string
                                Spectrum spec1 = sl.spectrum(spectrumIndex, true);
                                MZIntensityPairList peaks = new MZIntensityPairList();
                                spec1.getMZIntensityPairs(ref peaks);
                                Set<Peak> peakList = new Set<Peak>();
                                    
                                //get base peak and TIC
                                double basePeakValue = 0;
                                double TICValue = 0;
                                CVParamList list = spec1.cvParams;
                                foreach (CVParam CVP in list)
                                {
                                    if (CVP.name == "base peak intensity")
                                    {
                                        basePeakValue = CVP.value;
                                    }
                                    if (CVP.name == "total ion current")
                                    {
                                        TICValue = CVP.value;
                                    }
                                }
                                string basePeak = basePeakValue.ToString();
                                string TIC = TICValue.ToString();

                                //very important. Surendra put them here
                                //to change those with modifications {} into a format
                                //that fragment method will accept. 
                                string interpretation = vi.ToSimpleString();
                                    
                                Peptide peptide = new Peptide(interpretation, ModificationParsing.ModificationParsing_Auto, ModificationDelimiter.ModificationDelimiter_Brackets);
                                    
                                Fragmentation fragmentation = peptide.fragmentation(true, true);
                                //prepare the qualified peaklist
                                double intenThreshhold = 0;
                                int totalIntenClass = 0;

                                double[] intenArray = new double[peaks.Count];
                                //used during the foreach loop
                                int indexPeaks = 0;

                                //get all the peaks no matter how small they are
                                //then calculate the threshhold TIC
                                //test here
                                foreach (MZIntensityPair mzIntensity in peaks)
                                {
                                    //Peak p = new Peak(mzIntensity.mz, mzIntensity.intensity);
                                    //peakList.Add(p);
                                    intenArray[indexPeaks] = mzIntensity.intensity;
                                    indexPeaks++;
                                }
                                Array.Sort(intenArray);
                                Array.Reverse(intenArray, 0, peaks.Count);

                                //if currTIC>=cutoff, then break
                                double currTIC = 0;
                                double cutOffTIC = TicCutoffPercentage * TICValue;
                                foreach (double inten in intenArray)
                                {
                                    currTIC = currTIC + inten;
                                    if (currTIC < cutOffTIC)
                                    {
                                        intenThreshhold = inten;
                                        totalIntenClass++;
                                    }
                                    else break;
                                }


                                //then based on that, generate a new peaklist that contains only ABC peaks
                                //then calculate the intensity classes
                                foreach (MZIntensityPair mzIntensity in peaks)
                                {
                                    if (mzIntensity.intensity >= intenThreshhold)
                                    {
                                        //note 0 here. This is to tell people that the orbiorbi fragment charge is unknown. 
                                        //the peaklist will be updated later. 
                                        Peak p = new Peak(mzIntensity.mz, mzIntensity.intensity,0);
                                        peakList.Add(p);
                                    }
                                }
                                Console.WriteLine("nativeID: =============" + idOrIndex);
                                //update peaklist for charge states.
                                peakList = Package.chargeAssignment(peakList);
                                        //int ones = 0;
                                        //int twos = 0;
                                        //foreach (var peak in peakList)
                                        //{
                                        //    //Console.WriteLine(peak.fragmentCharge);
                                        //    if (peak.fragmentCharge == 1) ones++;
                                        //    else if (peak.fragmentCharge == 2) twos++;
                                        //}

                                        //Console.WriteLine("charge 1: " + ones);
                                        //Console.WriteLine("charge 2: " + twos);

                                //rowDic contains row information of each peptide bond
                                Dictionary<int, string> rowDic = new Dictionary<int, string>();
                                //intensityDic contains intensity information of fragment ions for each peptide bond
                                Dictionary<int, List<double>> intensityDic = new Dictionary<int, List<double>>();
                                //commonList contains the common intensities
                                List<double> duplicateList = new List<double>(); 

                                //call the method
                                List<double> completeIntensityList = new List<double>();
                                for (int k = 1; k < len; k++)
                                {
                                    List<double> intensityList = new List<double>();
                                    
                                    string bion = pepSequence.Substring(0, k);
                                    string yion = pepSequence.Substring(k, len - k);

                                    int NR = Package.parseAAResidues(bion, 'R');
                                    int NK = Package.parseAAResidues(bion, 'K');
                                    int NH = Package.parseAAResidues(bion, 'H');
                                    int NL = k;
                                    int CR = Package.parseAAResidues(yion, 'R');
                                    int CK = Package.parseAAResidues(yion, 'K');
                                    int CH = Package.parseAAResidues(yion, 'H');
                                    int CL = len - k;
                                    int pepBond = k;
                                    int NBasicAA = NR + NK + NH;
                                    int CBasicAA = CR + CK + CH;
                                    string AA = NR + "," + NK + "," + NH + "," + NL + "," + CR + "," + CK + "," + CH + "," + CL;
                                    
                                    double[] bCharge = new double[z + 1];
                                    double[] yCharge = new double[z + 1];
                                    //double[] bIntensity = new double[z + 1];
                                    //double[] yIntensity = new double[z + 1];
                                    //add variable for "real" charges
                                    
                                    if (model == 0) //naive model
                                    {
                                        int[] bFragmentCharge = new int[z + 1];
                                        int[] yFragmentCharge = new int[z + 1];
                                        //return the b ion charge 1 Intensity, if matched
                                        for (int i = 1; i < z; i++)
                                        {
                                            bCharge[i] = fragmentation.b(k, i);
                                            yCharge[i] = fragmentation.y(len - k, i);
                                            //change for ORBI-ORBI purposes.
                                            Peak bmatched = Package.findClose(peakList, bCharge[i], bCharge[i] * 30 * Math.Pow(10, -6));
                                            Peak ymatched = Package.findClose(peakList, yCharge[i], yCharge[i] * 30 * Math.Pow(10, -6));
                                            if (bmatched != null)
                                            {
                                                //bIntensity[i] = bmatched.rankOrIntensity;
                                                int fragmentCharge = bmatched.fragmentCharge;
                                                if (fragmentCharge == i)
                                                    bFragmentCharge[i] = 3;
                                                else bFragmentCharge[i] = 2;
                                                //intensityList.Add(bmatched.rankOrIntensity);
                                                //completeIntensityList.Add(bmatched.rankOrIntensity);
                                            }
                                            else bFragmentCharge[i] = 1;
                                            //else bIntensity[i] = 0;
                                            if (ymatched != null)
                                            {
                                                //yIntensity[i] = ymatched.rankOrIntensity;
                                                //yFragmentCharge[i] = ymatched.fragmentCharge;
                                                int fragmentCharge = ymatched.fragmentCharge;
                                                if (fragmentCharge == i)
                                                    yFragmentCharge[i] = 3;
                                                else yFragmentCharge[i] = 2;
                                                //intensityList.Add(ymatched.rankOrIntensity);
                                                //completeIntensityList.Add(ymatched.rankOrIntensity);
                                            }
                                            else yFragmentCharge[i] = 1;
                                            //else yIntensity[i] = 0;
                                        }
                                        string finalString = idOrIndex + "," + pepSequence + "," + pepBond + "," + bFragmentCharge[1] + "," + bFragmentCharge[2] + "," + yFragmentCharge[1] + "," + yFragmentCharge[2];
                                        file.WriteLine(finalString);
                                    }
                                    else if (model == 1) //my binary logistic regression basophile model
                                    {
                                        int b1=0, b2=0, y1=0, y2 = 0;
                                        double y1_logit = 0.1098112 * NR + 0.2085831 * NK + 0.1512109 * NH + 0.0460839 * NL
                                                        - 0.3872417 * CR - 0.3684911 * CK - 0.1634741 * CH - 0.1693931 * CL + 1.2632997;
                                        double y2_logit =-0.6345364 * NR - 0.3365917 * NK - 0.4577882 * NH - 0.1492703 * NL
                                                        + 0.7738133 * CR + 0.6036758 * CK + 0.5942542 * CH + 0.0701467 * CL + 0.0806280;
                                        double b1_logit = 0.0801432 * NR - 0.1088081 * NK - 0.1338220 * NH - 0.1413059 * NL
                                                        - 0.3157957 * CR - 0.2708274 * CK - 0.3703136 * CH + 0.0157418 * CL + 1.2124699;
                                        double b2_logit = 0.8606449 * NR + 0.2763119 * NK + 0.4969152 * NH + 0.0685712 * NL
                                                        - 1.3346995 * CR - 1.0977316 * CK - 1.0973677 * CH - 0.2028884 * CL + 1.9355980;
                                        if (b1_logit > -0.5)
                                        {
                                            double mz_b1 = fragmentation.b(k, 1);
                                            Peak matched = Package.findClose(peakList, mz_b1, mz_b1 * 30 * Math.Pow(10, -6));
                                            if (matched != null)
                                            {
                                                if (matched.fragmentCharge == 1) b1 = 3;
                                                else b1 = 2;
                                            }
                                            else b1 = 1;
                                        }
                                        else b1 = 0;
                                        if (b2_logit > 0)
                                        {
                                            double mz_b2 = fragmentation.b(k, 2);
                                            Peak matched = Package.findClose(peakList, mz_b2, mz_b2 * 30 * Math.Pow(10, -6));
                                            if (matched != null)
                                            {
                                                if (matched.fragmentCharge == 2) b2 = 3;
                                                else b2 = 2;
                                            }
                                            else b2 = 1;
                                        }
                                        else b2 = 0;
                                        if (y1_logit > -0.5)
                                        {
                                            double mz_y1 = fragmentation.y(len - k, 1);
                                            Peak matched = Package.findClose(peakList, mz_y1, mz_y1 * 30 * Math.Pow(10, -6));
                                            if (matched != null)
                                            {
                                                if (matched.fragmentCharge == 1) y1 = 3;
                                                else y1 = 2;
                                            }
                                            else y1 = 1;
                                        }
                                        else y1 = 0;
                                        if (y2_logit > -0.5)
                                        {
                                            double mz_y2 = fragmentation.y(len - k, 2);
                                            Peak matched = Package.findClose(peakList, mz_y2, mz_y2 * 30 * Math.Pow(10, -6));
                                            if (matched != null)
                                            {
                                                if (matched.fragmentCharge == 2) y2 = 3;
                                                else y2 = 2;
                                            }
                                            else y2 = 1;
                                        }
                                        else y2 = 0;
                                        string finalString = idOrIndex + "," + pepSequence + "," + pepBond + "," + b1 + "," + b2 + "," + y1 + "," + y2;
                                        file.WriteLine(finalString);
                                    }

                                    else if (model == 2) //Surendra's ordinal model
                                    {
                                        int b1 = 0, b2 = 0, y1 = 0, y2 = 0;
                                        double logit = NR * 0.9862 + NH * 0.8772 + NK * 0.7064 + NL * 0.4133 
                                                     - CR * 1.1688 - CH * 0.3948 - CK * 0.6710 - CL * 0.4859;
                                        //charge Label = "1", generate b+,y++
                                        if (logit < -2.2502)
                                        {
                                            double mz_b1 = fragmentation.b(k, 1);
                                            Peak b1matched = Package.findClose(peakList, mz_b1, mz_b1 * 30 * Math.Pow(10, -6));
                                            if (b1matched != null)
                                            {
                                                if (b1matched.fragmentCharge == 1) b1 = 3;  //matched, charge agree
                                                else b1 = 2; //matched, charge does not agree
                                            }
                                            else b1 = 1; //unmatched, but predicted. 

                                            double mz_y2 = fragmentation.y(len - k, 2);
                                            Peak y2matched = Package.findClose(peakList, mz_y2, mz_y2 * 30 * Math.Pow(10, -6));
                                            if (y2matched != null)
                                            {
                                                if (y2matched.fragmentCharge == 2) y2 = 3;  //matched, charge agree
                                                else y2 = 2; //matched, charge does not agree
                                            }
                                            else y2 = 1; //unmatched, but predicted. 
                                        }
                                        //ambiLabel = "2", generate b+, y+, b++, y++
                                        else if (logit < 0.7872)
                                        {
                                            double mz_b1 = fragmentation.b(k, 1);
                                            Peak b1matched = Package.findClose(peakList, mz_b1, mz_b1 * 30 * Math.Pow(10, -6));
                                            if (b1matched != null)
                                            {
                                                if (b1matched.fragmentCharge == 1) b1 = 3;  //matched, charge agree
                                                else b1 = 2; //matched, charge does not agree
                                            }
                                            else b1 = 1; //unmatched, but predicted. 

                                            double mz_y1 = fragmentation.y(len - k, 1);
                                            Peak y1matched = Package.findClose(peakList, mz_y1, mz_y1 * 30 * Math.Pow(10, -6));
                                            if (y1matched != null)
                                            {
                                                if (y1matched.fragmentCharge == 1) y1 = 3;  //matched, charge agree
                                                else y1 = 2; //matched, charge does not agree
                                            }
                                            else y1 = 1; //unmatched, but predicted.

                                            double mz_b2 = fragmentation.b(k, 2);
                                            Peak b2matched = Package.findClose(peakList, mz_b2, mz_b2 * 30 * Math.Pow(10, -6));
                                            if (b2matched != null)
                                            {
                                                if (b2matched.fragmentCharge == 2) b2 = 3;  //matched, charge agree
                                                else b2 = 2; //matched, charge does not agree
                                            }
                                            else b2 = 1; //unmatched, but predicted. 

                                            double mz_y2 = fragmentation.y(len - k, 2);
                                            Peak y2matched = Package.findClose(peakList, mz_y2, mz_y2 * 30 * Math.Pow(10, -6));
                                            if (y2matched != null)
                                            {
                                                if (y2matched.fragmentCharge == 2) y2 = 3;  //matched, charge agree
                                                else y2 = 2; //matched, charge does not agree
                                            }
                                            else y2 = 1; //unmatched, but predicted. 
                                        }
                                        //ambiLabel = "3", generate b++,y+
                                        else
                                        {
                                            double mz_b2 = fragmentation.b(k, 2);
                                            Peak b2matched = Package.findClose(peakList, mz_b2, mz_b2 * 30 * Math.Pow(10, -6));
                                            if (b2matched != null)
                                            {
                                                if (b2matched.fragmentCharge == 2) b2 = 3;  //matched, charge agree
                                                else b2 = 2; //matched, charge does not agree
                                            }
                                            else b2 = 1; //unmatched, but predicted. 

                                            double mz_y1 = fragmentation.y(len - k, 1);
                                            Peak y1matched = Package.findClose(peakList, mz_y1, mz_y1 * 30 * Math.Pow(10, -6));
                                            if (y1matched != null)
                                            {
                                                if (y1matched.fragmentCharge == 1) y1 = 3;  //matched, charge agree
                                                else y1 = 2; //matched, charge does not agree
                                            }
                                            else y1 = 1; //unmatched, but predicted. 
                                        }
                                        string finalString = idOrIndex + "," + pepSequence + "," + pepBond + "," + b1 + "," + b2 + "," + y1 + "," + y2;
                                        file.WriteLine(finalString);
                                    }
                                    
                                    //intensityDic.Add(pepBond, intensityList);

                                    
                                    //rowDic.Add(pepBond, finalString);
                                }//end for each peptide bond
                            }//end if z==3
                        }//end foreach peptide
            }//end using
        }//end idpReader

        public static void Main(string[] args)
        {
            string idpxml = "Z:\\home\\dwang\\fragmentation\\20091202-Scripps-Autophagy-OrbiOrbi\\human_baso_S\\113009_ATG7_WT_orbi-orbi.idpXML";
            string mzxml = "Z:\\home\\dwang\\fragmentation\\20091202-Scripps-Autophagy-OrbiOrbi\\113009_ATG7_WT_orbi-orbi.mzXML";
            double cutoff = 0.98;
            int z = 3;
            int model = 2; //0: naive; 1: logistic baso; 2: ordinal baso
            idpReader(idpxml, mzxml, cutoff, z, model);






            //concern: movenext in findclose function will mess up if not ordered correctly
            //Peak p1 = new Peak(3.000, 100.338);
            //Peak p2 = new Peak(4.007, 234390);
            //Peak p3 = new Peak(5.4, 3423.7834);
            //Peak p4 = new Peak(3.0, 4.0);
            //Peak p5 = new Peak(6.5, 402003.0923);
            //Peak p6 = new Peak(5.014, 100);
            //Set<Peak> peaklist = new Set<Peak>();
            //peaklist.Add(p1);
            //peaklist.Add(p2);
            //peaklist.Add(p3);
            //peaklist.Add(p4);
            //peaklist.Add(p5);
            //peaklist.Add(p6);
            //Console.WriteLine("starting......\n");
            //Set<Peak> temp = Package.chargeAssignment(peaklist);
            //foreach (var peak in temp)
            //{
            //    Console.WriteLine("mz= {0}; intensity={1}; z={2}\n", peak.mz, peak.rankOrIntensity, peak.fragmentCharge);
            //}
            
            
            // Keep the console window open in debug mode.
            Console.WriteLine("Press any key to exit.");
            Console.ReadKey();
        }

    }//end IdpReader
}
