using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using pwiz.CLI.analysis;
using pwiz.CLI.data;
using pwiz.CLI.msdata;
using pwiz.CLI.proteome;
using SourceList = System.Collections.Generic.Set<IDPicker.SourceInfo>;



namespace IDPReader
{
    class Package
    {
        IDPicker.Workspace workspace;

        public static void readWorkspace(string idpXMLPath)
        {
            try
            {
                Package pp = new Package();
                pp.workspace = new IDPicker.Workspace();
                StreamReader reader = new StreamReader(idpXMLPath);
                pp.workspace.readPeptidesXml(reader, "", 0.05f, 1);
            }
            catch (Exception e)
            {
                Console.WriteLine("Error while parsing idpXML");
                Console.Write(e.StackTrace);
            }
        }

        //how to assemble the workspace from an idpXML
        public static void loadWorkspace(ref IDPicker.Workspace workspace, string assemblyFile)
        {
            if (workspace == null)
            {
                return;
            }

            try
            {
                workspace.readPeptidesXml(new StreamReader(assemblyFile), "", 1.0f, 1);
            }
            catch (Exception e)
            {
                Console.Error.WriteLine(e.StackTrace.ToString());
                Console.Error.WriteLine("Error reading input filepath \"" + assemblyFile + "\": " + e.Message);
                Environment.Exit(1);
            }

            // Create peptide groups by grouping peptide identifications
            // that mapped to same proteins into a single group.

            workspace.assemblePeptideGroups();

        }

        /// <summary>
        ///  This function takes a set of peaks, a lookup mass, and 
        ///  finds the closest peak with in a certain tolerance. If 
        ///  multiple peaks exists with in the window, most intense peak
        ///  is selected
        /// </summary>
        /// <param name="peaks">A spectrum</param>
        /// <param name="mz">Look up m/z</param>
        /// <param name="tolerance">Mass tolerance for look-up</param>
        /// <returns>Peak if found</returns>
        public static Peak findNear(Set<Peak> peaks, double mz, double tolerance)
        {
            Set<Peak>.Enumerator cur, min, max;

            min = peaks.LowerBound(new Peak(mz - tolerance));
            max = peaks.LowerBound(new Peak(mz + tolerance));
            if (!min.IsValid && !max.IsValid)
                return null;
            if (!min.IsValid && max.IsValid)
                return max.Current;
            if (min.IsValid && !max.IsValid)
                return min.Current;
            if (min.Current == max.Current)
                return null;

            // If we found multiple matching peaks, 
            // return the peak with best intensity.
            Peak best = min.Current;
            double bestIntensityOrRank = best.rankOrIntensity;
            for (cur = min; cur.Current != max.Current; cur.MoveNext())
            {
                double curRank = cur.Current.rankOrIntensity;
                if (curRank > bestIntensityOrRank)
                {
                    bestIntensityOrRank = curRank;
                    best = cur.Current;
                }
            }
            return best;
        }

        /// <summary>
        ///  This function takes a set of peaks, a lookup mass, and 
        ///  finds the closest peak with in a certain tolerance. If 
        ///  multiple peaks exists with in the window, most closest peak
        ///  is selected
        /// </summary>
        /// <param name="peaks">A spectrum</param>
        /// <param name="mz">Look up m/z</param>
        /// <param name="tolerance">Mass tolerance for look-up</param>
        /// <returns>Peak if found</returns>
        public static Peak findClose(Set<Peak> peaks, double mz, double tolerance)
        {
            Set<Peak>.Enumerator cur, min, max;

            min = peaks.LowerBound(new Peak(mz - tolerance));
            max = peaks.LowerBound(new Peak(mz + tolerance));
            if (!min.IsValid && !max.IsValid)
                return null;
            if (!min.IsValid && max.IsValid)
                return max.Current;
            if (min.IsValid && !max.IsValid)
                return min.Current;
            if (min.Current == max.Current)
                return null;

            // If we found multiple matching peaks, 
            // return the closest peak.
            Peak best = min.Current;

            //find the peak closest to the desired mz
            double minDiff = Math.Abs(mz - best.mz);
            for (cur = min; cur.Current != max.Current; cur.MoveNext())
            {
                
                double curDiff = Math.Abs(mz - cur.Current.mz);
                if (curDiff < minDiff)
                {
                    minDiff = curDiff;
                    best = cur.Current;
                }
            }
            return best;
        }

        /// <summary>
        ///  This function takes a set of peaks, iterate through the peaks 
        ///  and assign each fragment with a charge label with:
        ///  0: unknown; 1: +1; 2: +2; 3:+3 et al. 
        ///  note: the spectrum peaks comes from a +3 precursor. 
        ///  also I require 5% basepeak intensity cutoff for the peaklist
        ///  BASED ON FINDCLOSE METHOD. 
        /// </summary>
        public static Set<Peak> chargeAssignment(Set<Peak> peaklist)
        {
            Set<Peak> finalPeakList = new Set<Peak>();
            //first step should be sort the peaks from least to highest by mz
            Set<Peak>.Enumerator cur, min, max, prev;
            //find the min mz and max mz
            double minMZ = 0;
            //double intensity_min = 0;
            double maxMZ = 0;
            //double intensity_max = 0;
            foreach (var peak in peaklist)
            {
                if (peak.mz > maxMZ)
                {
                    maxMZ = peak.mz;
                    //intensity_max = peak.rankOrIntensity;
                }
                if (peak.mz < minMZ)
                {
                    minMZ = peak.mz;
                    //intensity_min = peak.rankOrIntensity;
                }
            }
            min = peaklist.LowerBound(new Peak(minMZ));
            max = peaklist.LowerBound(new Peak(maxMZ));

            //find the difference of adjacent peaks
            int size = peaklist.Count;
            double prev_mz = min.Current.mz;
            //double prev_diff = 0;

            int[] assignment = new int[size - 1];
            int itr = 0;
            for (cur = min; cur.Current != max.Current; cur.MoveNext())
            {
                //the diff of the first element (min) is 0
                //Console.WriteLine("current mz is: " + cur.Current.mz);
                //double 1/2H = 0.503638; round to 0.504
                //double H = 1.007276; round to 1.007
                //it turns out that TWO peaks will be enough for 
                double diff = Math.Round(cur.Current.mz - prev_mz,1);
                //Console.WriteLine(diff);
                if (diff == 0.5)
                {
                    //assignment[itr - 2] = 2;
                    assignment[itr - 1] = 2;
                    assignment[itr] = 2;
                }
                else if (diff == 1)
                {
                    //assignment[itr - 2] = 1;
                    assignment[itr - 1] = 1;
                    assignment[itr] = 1;
                }
                else assignment[itr] = 0;
                //Console.WriteLine("========\n");
                //Console.WriteLine(prev_diff + "\n");
                //Console.WriteLine(diff + "\n");
                prev = cur;
                //prev_diff = diff;
                prev_mz = cur.Current.mz;
                itr++;
            }

            //if the diff patterns shows a successive 3-peak of 0.5 (which is 2 diff), then assign 2
            //if the diff goes like 0.33333, then assign 3 or does not count?
            //if the diff goes like 1, then assign 1. 
            //otherwise assign just 0
            //the key is 5% cutoff. 
            
            min = peaklist.LowerBound(new Peak(minMZ));
            max = peaklist.LowerBound(new Peak(maxMZ));
            //Console.WriteLine(min.Current.mz);
            //Console.WriteLine(max.Current.mz);
            itr = 0;
            for (cur = min; cur.Current != max.Current; cur.MoveNext())
            {
                //Console.WriteLine("adding peak\n");
                finalPeakList.Add(new Peak(cur.Current.mz, cur.Current.rankOrIntensity, assignment[itr]));
                //Console.WriteLine(assignment[itr]);
                itr++;
            }
            return finalPeakList;
        }



        //given a fragment ion sequence, return the amino acid composition (number of aa)
        public static int parseAAResidues(string pepSequence, char AA)
        {
            int num = 0;

            char[] charArr = pepSequence.ToCharArray();

            foreach (char c in charArr)
            {
                if (c == AA) num++;

            }
            return num;
        }

        //find replicate
        public static List<double> findCommon(List<double> inputList)
        {
            Dictionary<double, int> uniqueStore = new Dictionary<double, int>();
            List<double> finalList = new List<double>();

            foreach (double currValue in inputList)
            {
                if (!uniqueStore.ContainsKey(currValue))
                {
                    uniqueStore.Add(currValue, 0);
                }
                else finalList.Add(currValue);
            }
            return finalList;
        } 


        ///<summary>
        ///original code for grabbing the required ions
        ///get pepSequence,source,scanID, write into a csv file. 
        ///have tons of information: theretical ion intensity, pep info, chargeLabel...
        ///</summary>
        
        public static List<string> idpReader_original(string idpXMLFile, string mzMLFile, double TicCutoffPercentage, int z, List<string> pepList, List<string> output)
        {
            //get the path and filename of output csv file:
            string fileName = Path.GetFileNameWithoutExtension(idpXMLFile);
            string filePath = Path.GetDirectoryName(idpXMLFile);
            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, idpXMLFile);


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

                            //make sure that the peptide is what we want
                            if (pepList.Contains(pepSequence))
                            {
                                // Look up the index with nativeID
                                object idOrIndex = null;
                                if (sItr.Value.nativeID != null && sItr.Value.nativeID.Length > 0)
                                    idOrIndex = sItr.Value.nativeID;
                                int spectrumIndex = sl.find(idOrIndex as string);
                                // Trust the local index, if the nativeID lookup fails
                                if (spectrumIndex >= sl.size())
                                    spectrumIndex = sItr.Value.id.index;
                                // Bail of the loca index is larger than the spectrum list size
                                if (spectrumIndex >= sl.size())
                                    throw new Exception("Can't find spectrum associated with the index.");

                                //Console.WriteLine(idOrIndex.ToString());
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
                                        Peak p = new Peak(mzIntensity.mz, mzIntensity.intensity);
                                        peakList.Add(p);
                                    }
                                }


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
                                    int R = NR - CR;
                                    int K = NK - CK;
                                    int H = NH - CH;
                                    int L = NL - CL;
                                    int pepBond = k;
                                    int NBasicAA = NR + NK + NH;
                                    int CBasicAA = CR + CK + CH;
                                    string AA = NBasicAA + "," + CBasicAA + "," + NR + "," + NK + "," + NH + "," + NL + "," + CR + "," + CK + "," + CH + "," + CL + "," + R + "," + K + "," + H + "," + L;


                                    double[] bCharge = new double[z + 1];
                                    double[] yCharge = new double[z + 1];
                                    double[] bIntensity = new double[z + 1];
                                    double[] yIntensity = new double[z + 1];
                                    string bIonIntensity = "";
                                    string yIonIntensity = "";
                                    //to judge if the sum of intensities are 0
                                    //so to exclude the case with all "0s"
                                    double sumIntensity = 0;

                                    //return the b ion charge 1 Intensity, if matched
                                    for (int i = 1; i <= z; i++)
                                    {
                                        bCharge[i] = fragmentation.b(k, i);
                                        yCharge[i] = fragmentation.y(len - k, i);
                                        //change for Q-star purposes.
                                        Peak bmatched = Package.findClose(peakList, bCharge[i], 70 * bCharge[i] * Math.Pow(10, -6));
                                        Peak ymatched = Package.findClose(peakList, yCharge[i], 70 * yCharge[i] * Math.Pow(10, -6));
                                        if (bmatched != null)
                                        {
                                            bIntensity[i] = bmatched.rankOrIntensity;
                                            intensityList.Add(bmatched.rankOrIntensity);
                                            completeIntensityList.Add(bmatched.rankOrIntensity);
                                        }
                                        else bIntensity[i] = 0;
                                        if (ymatched != null)
                                        {
                                            yIntensity[i] = ymatched.rankOrIntensity;
                                            intensityList.Add(ymatched.rankOrIntensity);
                                            completeIntensityList.Add(ymatched.rankOrIntensity);
                                        }
                                        else yIntensity[i] = 0;

                                        sumIntensity = sumIntensity + bIntensity[i] + yIntensity[i];

                                        //record b/y ion intensity information into a string
                                        bIonIntensity = bIonIntensity + "," + bIntensity[i];
                                        yIonIntensity = yIonIntensity + "," + yIntensity[i];

                                    }

                                    intensityDic.Add(pepBond, intensityList);
                                    //to determine charge label, need to split by precursor charge
                                    //first need to make a metric to determine if all intensities are "0"

                                    if (z == 3)
                                    {
                                        if (sumIntensity != 0)
                                        {
                                            ////////////////////////////////////////////////
                                            //set the ambiguity label as follows:
                                            //-3: (0/+3) y3 only
                                            //-2: (0/+3, +1/+2) y3, b1y2
                                            //-1: (+1/+2): b1y2
                                            //0: (+1/+2, +2/+1): b1y2, b2y1
                                            //1: (+2/+1): b2y1
                                            //2: (+2/+1, +3/0): b2y1, b3
                                            //3: (+3/0): b3 only
                                            ////////////////////////////////////////////////

                                            double b1 = bIntensity[1];
                                            double b2 = bIntensity[2];
                                            double y1 = yIntensity[1];
                                            double y2 = yIntensity[2];
                                            double b3 = bIntensity[3];
                                            double y3 = yIntensity[3];
                                            double b1y2 = b1 + y2;
                                            double b2y1 = b2 + y1;
                                            string ambiguityLabel = "";
                                            //first part: set the intensity group: y3, b1y2, b2y1, b3
                                            //if one group was found, set the label
                                            //if two were found, but adjacent to each other, then ambiguity label is set
                                            if (y3 != 0 && (b1y2 + b2y1 + b3) == 0) ambiguityLabel = "-3";
                                            else if (y3 != 0 && b1y2 != 0 && (b2y1 + b3) == 0) ambiguityLabel = "-2";
                                            else if (b1y2 != 0 && (y3 + b2y1 + b3) == 0) ambiguityLabel = "-1";
                                            else if (b1y2 != 0 && b2y1 != 0 && (y3 + b3) == 0) ambiguityLabel = "0";
                                            else if (b2y1 != 0 && (y3 + b1y2 + b3) == 0) ambiguityLabel = "1";
                                            else if (b2y1 != 0 && b3 != 0 && (y3 + b1y2) == 0) ambiguityLabel = "2";
                                            else if (b3 != 0 && (y3 + b1y2 + b2y1) == 0) ambiguityLabel = "3";
                                            else ambiguityLabel = "error";

                                            string finalString = idOrIndex + "," + pepSequence + "," + basePeak + "," + TIC + bIonIntensity + yIonIntensity + "," + len + "," + pepBond + "," + AA + "," + ambiguityLabel;
                                            rowDic.Add(pepBond, finalString);
                                        }
                                    }
                                }//end for each peptide bond

                                //now we have: rowDic, intensityDic for each pep bond
                                //and we have: and completeIntensityList for each peptide
                                //the purpose of this is to remove such rows with duplicate matches
                                duplicateList = Package.findCommon(completeIntensityList);
                                foreach (int bond in rowDic.Keys)
                                {
                                    bool unique = true;
                                    foreach (double inten in duplicateList)
                                    {
                                        if (intensityDic[bond].Contains(inten))
                                        {
                                            unique = false;
                                            Console.WriteLine("kick");
                                            break;
                                        }
                                    }
                                    if (unique)
                                    {
                                        output.Add(rowDic[bond]);
                                    }
                                }
                            }//end of if peplist contains pepsequence
                        }//end if z==3
                    }//end foreach peptide
            return output;
        }//end idpReader

        
        ///<summary>
        ///secondary version grabbing the required ions
        ///what is the difference of original and this version?
        ///this version does grab fragments, yet only 5% cutoff
        ///does not mimic what myrimatch does
        /// </summary>
        public static List<string> idpReader(string idpXMLFile, string mzMLFile, double TicCutoffPercentage, int z)     
        {
            //test
            //Console.WriteLine("test1");
            int repeatitive = 0;

            List<string> output = new List<string>();
            //get the path and filename of output csv file:
            string fileName = Path.GetFileNameWithoutExtension(idpXMLFile);
            string filePath = Path.GetDirectoryName(idpXMLFile);
            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, idpXMLFile);

           // Console.WriteLine("test2");

            MSDataFile foo = new MSDataFile(mzMLFile);
            SpectrumList sl = foo.run.spectrumList;

            //Console.WriteLine("test3");
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

                            //make sure that the peptide is what we want
                            
                                // Look up the index with nativeID
                                object idOrIndex = null;
                                if (sItr.Value.nativeID != null && sItr.Value.nativeID.Length > 0)
                                    idOrIndex = sItr.Value.nativeID;
                                int spectrumIndex = sl.find(idOrIndex as string);
                                // Trust the local index, if the nativeID lookup fails
                                if (spectrumIndex >= sl.size())
                                    spectrumIndex = sItr.Value.id.index;

                                //test
                                //Console.WriteLine(spectrumIndex);
                                //// Bail of the loca index is larger than the spectrum list size
                                //if (spectrumIndex >= sl.size())
                                //    throw new Exception("Can't find spectrum associated with the index.");

                                //Console.WriteLine(idOrIndex.ToString());
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
                                double basePeakCutOff = 0.05 * basePeakValue;
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
                                        Peak p = new Peak(mzIntensity.mz, mzIntensity.intensity);
                                        peakList.Add(p);
                                    }
                                }


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
                                    int R = NR - CR;
                                    int K = NK - CK;
                                    int H = NH - CH;
                                    int L = NL - CL;
                                    int pepBond = k;
                                    int NBasicAA = NR + NK + NH;
                                    int CBasicAA = CR + CK + CH;
                                    string AA = NR + "," + NK + "," + NH + "," + NL + "," + CR + "," + CK + "," + CH + "," + CL;


                                    double[] bCharge = new double[z + 1];
                                    double[] yCharge = new double[z + 1];
                                    double[] bIntensity = new double[z + 1];
                                    double[] yIntensity = new double[z + 1];
                                    string bIonIntensity = "";
                                    string yIonIntensity = "";
                                    //to judge if the sum of intensities are 0
                                    //so to exclude the case with all "0s"
                                    double sumIntensity = 0;

                                    //return the b ion charge 1 Intensity, if matched
                                    for (int i = 1; i <= z; i++)
                                    {
                                        bCharge[i] = fragmentation.b(k, i);
                                        yCharge[i] = fragmentation.y(len - k, i);
                                        //change for LTQ/ORBI.
                                        Peak bmatched = Package.findClose(peakList, bCharge[i], 0.5);
                                        Peak ymatched = Package.findClose(peakList, yCharge[i], 0.5);
                                        if (bmatched != null && bmatched.rankOrIntensity > basePeakCutOff)
                                        {
                                            bIntensity[i] = bmatched.rankOrIntensity;
                                            intensityList.Add(bmatched.rankOrIntensity);
                                            completeIntensityList.Add(bmatched.rankOrIntensity);
                                        }
                                        else bIntensity[i] = 0;
                                        if (ymatched != null && ymatched.rankOrIntensity > basePeakCutOff)
                                        {
                                            yIntensity[i] = ymatched.rankOrIntensity;
                                            intensityList.Add(ymatched.rankOrIntensity);
                                            completeIntensityList.Add(ymatched.rankOrIntensity);
                                        }
                                        else yIntensity[i] = 0;

                                        sumIntensity = sumIntensity + bIntensity[i] + yIntensity[i];

                                        //record b/y ion intensity information into a string
                                        bIonIntensity = bIonIntensity + "," + bIntensity[i];
                                        yIonIntensity = yIonIntensity + "," + yIntensity[i];

                                    }

                                    intensityDic.Add(pepBond, intensityList);
                                    //to determine charge label, need to split by precursor charge
                                    //first need to make a metric to determine if all intensities are "0"

                                    if (z == 2)
                                    {
                                        if (sumIntensity != 0)
                                        {
                                            double b1 = bIntensity[1];
                                            double b2 = bIntensity[2];
                                            double y1 = yIntensity[1];
                                            double y2 = yIntensity[2];
                                            string b1label = "0";
                                            string y1label = "0";

                                            if (b1 != 0) b1label = "1";
                                            if (y1 != 0) y1label = "1";

                                            string labels = b1label + "," + y1label;

                                            string finalString = AA + "," + labels;
                                            rowDic.Add(pepBond, finalString);
                                        }
                                    }

                                    if (z == 3)
                                    {
                                        if (sumIntensity != 0)
                                        {
                                            double b1 = bIntensity[1];
                                            double b2 = bIntensity[2];
                                            double y1 = yIntensity[1];
                                            double y2 = yIntensity[2];
                                            double b3 = bIntensity[3];
                                            double y3 = yIntensity[3];
                                            string b1label = "0";
                                            string b2label = "0";
                                            string y1label = "0";
                                            string y2label = "0";

                                            if (b1 != 0) b1label = "1";
                                            if (b2 != 0) b2label = "1";
                                            if (y1 != 0) y1label = "1";
                                            if (y2 != 0) y2label = "1";

                                            string labels = b1label + "," + b2label + "," + y1label + "," + y2label;

                                            string finalString = AA + "," + labels;
                                            rowDic.Add(pepBond, finalString);
                                        }
                                    }


                                }//end for each peptide bond

                                //now we have: rowDic, intensityDic for each pep bond
                                //and we have: and completeIntensityList for each peptide
                                //the purpose of this is to remove such rows with duplicate matches
                                duplicateList = Package.findCommon(completeIntensityList);
                                
                                foreach (int bond in rowDic.Keys)
                                {
                                    bool unique = true;
                                    foreach (double inten in duplicateList)
                                    {
                                        if (intensityDic[bond].Contains(inten))
                                        {
                                            unique = false;
                                            repeatitive++;
                                            break;
                                        }
                                    }
                                    if (unique)
                                    {
                                        output.Add(rowDic[bond]);
                                    }
                                }
                        }//end of boolcharge
                    }//end foreach peptide
            Console.WriteLine(repeatitive);
            return output;
        }//end idpReader

    }
}
