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



namespace Application
{
    public class Peak : IComparable<Peak>
    {

        public double mz;
        public double rankOrIntensity;

        public Peak(double mass, double rnkOrIntens)
        {
            mz = mass;
            rankOrIntensity = rnkOrIntens;
        }

        public Peak(Peak pk)
        {
            mz = pk.mz;
            rankOrIntensity = pk.rankOrIntensity;
        }

        public Peak(double mass)
        {
            mz = mass;
            rankOrIntensity = -1;
        }

        public int CompareTo(Peak rhs)
        {
            if (mz < rhs.mz)
                return -1;
            else if (mz > rhs.mz)
                return 1;
            return 0;
        }
    }

    class IdpReader
    {
        IDPicker.Workspace workspace;

        public void readWorkspace(string idpXMLPath)
        {
            try
            {
                workspace = new IDPicker.Workspace();
                StreamReader reader = new StreamReader(idpXMLPath);
                workspace.readPeptidesXml(reader, "", 0.05f, 1);
            }
            catch (Exception e)
            {
                Console.WriteLine("Error while parsing idpXML");
                Console.Write(e.StackTrace);
            }
        }


        //how to assemble the workspace from an idpXML
        public void loadWorkspace(ref IDPicker.Workspace workspace, string assemblyFile)
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

       
        //get pepSequence,source,scanID, write into a csv file. 
        //have tons of information: theretical ion intensity, pep info, chargeLabel...
        public static void idpReader(string idpXMLFile, string mzMLFile, double TicCutoffPercentage, int z)
        {
            //get the path and filename of output csv file:
            string fileName = Path.GetFileNameWithoutExtension(idpXMLFile);
            string filePath = Path.GetDirectoryName(idpXMLFile);
            string csvFile = Path.Combine(filePath, fileName) + "_" + z.ToString() + ".csv";
            
            IDPicker.Workspace workspace = new IDPicker.Workspace();
            IdpReader pp = new IdpReader();
            pp.loadWorkspace(ref workspace, idpXMLFile);

            using (StreamWriter file = new StreamWriter(csvFile))
            {
                if (z == 2)
                {
                    file.WriteLine("nativeID,pepSequence,basePeak,TIC,b1,b2,y1,y2,len,bond,NBasic,CBasic,NR,NK,NH,NL,CR,CK,CH,CL,ambiguityLabel");
                }
                if (z == 3)
                {
                    file.WriteLine("nativeID,pepSequence,basePeak,TIC,b1,b2,b3,y1,y2,y3,len,bond,NBasic,CBasic,NR,NK,NH,NL,CR,CK,CH,CL,ambiguityLabel");
                }
                
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
                                
                                Spectrum spec1 = sl.spectrum(spectrumIndex, true);
                                MZIntensityPairList peaks = new MZIntensityPairList();
                                spec1.getMZIntensityPairs(ref peaks);
                                Set<Peak> peakList = new Set<Peak>();
                                //get base peak and TIC and converted to string

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
                                string interpretation = vi.ToSimpleString();
                                Peptide peptide = new Peptide(interpretation, ModificationParsing.ModificationParsing_Auto, ModificationDelimiter.ModificationDelimiter_Brackets); 
                                
                                Fragmentation fragmentation = peptide.fragmentation(true, true);
                                //prepare the qualified peaklist
                                Console.WriteLine("prepare peaklist ->...");
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

                                //call the method
                                Console.WriteLine("fragment ion ->...");
                                for (int k = 1; k < len; k++)
                                {
                                    string bion = pepSequence.Substring(0, k);
                                    string yion = pepSequence.Substring(k, len - k);

                                    int NR = parseAAResidues(bion, 'R');
                                    int NK = parseAAResidues(bion, 'K');
                                    int NH = parseAAResidues(bion, 'H');
                                    int NL = k;
                                    int CR = parseAAResidues(yion, 'R');
                                    int CK = parseAAResidues(yion, 'K');
                                    int CH = parseAAResidues(yion, 'H');
                                    int CL = len - k;
                                    int R = NR - CR;
                                    int K = NK - CK;
                                    int H = NH - CH;
                                    int L = NL - CL;
                                    double pepBond = k;
                                    int NBasicAA = NR + NK + NH;
                                    int CBasicAA = CR + CK + CH;
                                    string AA = NBasicAA.ToString() + "," + CBasicAA.ToString() + "," + NR.ToString() + "," + NK.ToString() + "," + NH.ToString() + "," + NL.ToString() + "," + CR.ToString() + "," + CK.ToString() + "," + CH.ToString() + "," + CL.ToString();

                                    
                                    double[] bCharge = new double [z+1];
                                    double[] yCharge = new double [z+1];
                                    double[] bIntensity = new double [z+1];
                                    double[] yIntensity = new double [z+1];
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
                                        Peak bmatched = findNear(peakList, bCharge[i], 0.5);
                                        Peak ymatched = findNear(peakList, yCharge[i], 0.5);
                                        if (bmatched != null) bIntensity[i] = bmatched.rankOrIntensity;
                                        else bIntensity[i] = 0;
                                        if (ymatched != null) yIntensity[i] = ymatched.rankOrIntensity;
                                        else yIntensity[i] = 0;

                                        sumIntensity = sumIntensity + bIntensity[i] + yIntensity[i];

                                        //record b/y ion intensity information into a string
                                        bIonIntensity = bIonIntensity + "," + bIntensity[i];
                                        yIonIntensity = yIonIntensity + "," + yIntensity[i];

                                    }

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
                                            double b1y1 = b1 + y1;
                                            string ambiguityLabel = "";
                                            //string chargeLabel = "";
                                            //to judge the chargeLabel, there got be basic AA requirement as follows:
                                            //for chargeLabel -1 (0/+2), CBasicAA >=1
                                            //for chargeLabel +1 (+2/0), NBasicAA >=1
                                            //the following are not compatible with ambiguity label
                                            /*
                                            if (b1 > b2 && b1 > y2) chargeLabel = "0";
                                            else if (y1 > b2 && y1 > y2) chargeLabel = "0";
                                            else if (b2 > b1 && b2 > y1 && b2 > y2)
                                            {
                                                if (NBasicAA < 1) chargeLabel = "error";
                                                else chargeLabel = "1";
                                            }
                                            else if (y2 > b1 && y2 > b2 && y2 > y1)
                                            {
                                                if (CBasicAA < 1) chargeLabel = "error";
                                                else chargeLabel = "-1";
                                            }
                                            else chargeLabel = "empty";
                                             */

                                            //to judge ambiguity label
                                            //put b1+y1, b2, y3 for comparison, record the corresponding label
                                            //numbers: -2(0/+2), -1(0/+2, +1/+1), 0(+1/+1), 1(+1/+1, +2/0), 2(+2/0)
                                            if (b2 != 0 && (b1y1 + y2) == 0) 
                                            {
                                                if (NBasicAA < 1) ambiguityLabel = "FAA";
                                                else ambiguityLabel = "2";
                                            }
                                            else if (b2 != 0 && b1y1 != 0 && (y2 == 0)) 
                                            {
                                                if (NBasicAA < 1) ambiguityLabel = "faa";
                                                else ambiguityLabel = "1";
                                            }
                                            else if (b1y1 != 0 && (b2 + y2) == 0) ambiguityLabel = "0";
                                            else if (b1y1 != 0 && y2 != 0 && b2 == 0) 
                                            {
                                                if (CBasicAA < 1) ambiguityLabel = "faa";
                                                else ambiguityLabel = "-1";
                                            }
                                            else if (y2 != 0 && (b1y1 + b2) == 0) 
                                            {
                                                if (CBasicAA < 1) ambiguityLabel = "faa";
                                                else ambiguityLabel = "-2";
                                            }
                                            //here there are two cases excluded
                                            //b2+y2, b2 b1y1 y2. 
                                            else ambiguityLabel = "error";

                                            //here note that bIonIntensity string is like ",b1,b2,b3"...
                                            string finalString = idOrIndex + "," + pepSequence + "," + basePeak + "," + TIC + bIonIntensity + yIonIntensity + "," + len + "," + pepBond + "," + AA + "," +  ambiguityLabel;
                                            Console.WriteLine("seqIntensityList successfully added ->...");
                                            file.WriteLine(finalString);
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
                                            double b1y2 = b1 + y2;
                                            double b2y1 = b2 + y1;
                                            //string chargeLabel = "";
                                            string ambiguityLabel = "";

                                            //to judge the chargeLabel, there got be basic AA requirement as follows:
                                            //for chargeLabel -1 (+1/+2), CBasicAA >=1
                                            //for chargeLabel +1 (+2/+1), NBasicAA >=1
                                            //for chargeLabel -2 (+0/+3), CBasicAA >=2
                                            //for chargeLabel +2 (+3/+0), NBasicAA >=2
                                             //the following are not compatible with ambiguity label
                                            /*
                                            if (b1y2 > b2y1 && b1y2 > b3 && b1y2 > y3)
                                            {
                                                if (CBasicAA < 1)
                                                {
                                                    chargeLabel = "error";
                                                }
                                                else 
                                                {
                                                    chargeLabel = "-1";
                                                }
                                            }
                                            else if (b2y1 > b1y2 && b2y1 > b3 && b2y1 > y3)
                                            {
                                                if (NBasicAA < 1)
                                                {
                                                    chargeLabel = "error";
                                                }
                                                else
                                                {
                                                    chargeLabel = "1";
                                                }
                                            }
                                            else if (b3 > b1y2 && b3 > b2y1 && b3 > y3)
                                            {
                                                if (NBasicAA < 2)
                                                {
                                                    chargeLabel = "error";
                                                }
                                                else
                                                {
                                                    chargeLabel = "2";
                                                }
                                            }
                                            else if (y3 > b1y2 && y3 > b2y1 && y3 > b3)
                                            {
                                                if (CBasicAA < 2)
                                                {
                                                    chargeLabel = "error";
                                                }
                                                else
                                                {
                                                    chargeLabel = "-2";
                                                }
                                            }
                                            else
                                            {
                                                chargeLabel = "empty";
                                            }
                                             */


                                            //set the ambiguity label as follows:
                                            //-3: (0/+3) y3 only
                                            //-2: (0/+3, +1/+2) y3, b1y2
                                            //-1: (+1/+2): b1y2
                                            //0: (+1/+2, +2/+1): b1y2, b2y1
                                            //1: (+2/+1): b2y1
                                            //2: (+2/+1, +3/0): b2y1, b3
                                            //3: (+3/0): b3
                                            if (y3 != 0 && (b1y2 + b2y1 + b3) == 0) 
                                            {
                                                if (CBasicAA < 2) ambiguityLabel = "faa";
                                                else ambiguityLabel = "-3";
                                            }
                                            //if c-basic aa is more than 2, then it could be y3/y2
                                            //if c-basic aa is 1, then it could only be y2 (exclude y3)
                                            else if (y3 != 0 && b1y2 != 0 && (b2y1 + b3) == 0) 
                                            {
                                                if (CBasicAA >= 2) ambiguityLabel = "-2";
                                                else if (CBasicAA == 1) ambiguityLabel = "-1";
                                                else ambiguityLabel = "faa";
                                            }
                                            else if (b1y2 != 0 && (y3 + b2y1 + b3) == 0) 
                                            {
                                                if (CBasicAA < 1) ambiguityLabel = "faa";
                                                else ambiguityLabel = "-1";
                                            }
                                            else if (b1y2 != 0 && b2y1 != 0 && (y3 + b3) == 0) 
                                            {
                                                if (CBasicAA == 0 && NBasicAA > 0) ambiguityLabel = "1";
                                                else if (NBasicAA == 0 && CBasicAA > 0) ambiguityLabel = "-1";
                                                else if (NBasicAA > 0 && CBasicAA > 0) ambiguityLabel = "0";
                                                else ambiguityLabel = "faa";
                                            }
                                            else if (b2y1 != 0 && (y3 + b1y2 + b3) == 0) 
                                            {
                                                if (NBasicAA == 0) ambiguityLabel = "faa";
                                                else ambiguityLabel = "1";
                                            }
                                            else if (b2y1 != 0 && b3 != 0 && (y3 + b1y2) == 0) 
                                            {
                                                if (NBasicAA >= 2) ambiguityLabel = "2";
                                                else if (NBasicAA == 1) ambiguityLabel = "1";
                                                else ambiguityLabel = "faa";
                                            }
                                            else if (b3 != 0 && (y3 + b1y2 + b2y1) == 0) 
                                            {
                                                if (NBasicAA < 2) ambiguityLabel = "faa";
                                                else ambiguityLabel = "3";
                                            }
                                            else ambiguityLabel = "error";

                                            string finalString = idOrIndex + "," + pepSequence + "," + basePeak + "," + TIC + bIonIntensity + yIonIntensity + "," + len + "," + pepBond + "," + AA + "," + ambiguityLabel;
                                            Console.WriteLine("seqIntensityList successfully added ->...");
                                            file.WriteLine(finalString);
                                        }
                                        
                                    }
                 
                                }//end for
                            }//end if
                        }//end foreach
            }//end using
        }//end idpReader


        //the following is the flow chart of my current job
        ///////////////////////////////////////////////////
        //
        //     TASK FIRED
        //
        ///////////////////////////////////////////////////


        
        //
        //read idpxml get idlist
        //read pepxml, get all the ids which rank = 1
        //analyzing the class intensity class within the method, return a string
        //attach the string to the entire dictionary value
        public static Dictionary<string, string> pepXMLReader(string idpXMLFile, string pepXMLFile, string mzMLFile, string model, double TicCutoffPercentage, int z)
        {
            //first get IDs from idpXML file
            List<string> nativeIDList = new List<string>();
            IDPicker.Workspace workspace = new IDPicker.Workspace();
            IdpReader pp = new IdpReader();
            pp.loadWorkspace(ref workspace, idpXMLFile);


            MSDataFile foo = new MSDataFile(mzMLFile);
            SpectrumList sl = foo.run.spectrumList;

            foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
                foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
                    foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
                    {
                        IDPicker.ResultInstance ri = sItr.Value.results[1];
                        IDPicker.VariantInfo vi = ri.info.peptides.Min;

                        bool boolCharge = sItr.Value.id.charge.Equals(z);
                        if (boolCharge)
                        {
                            // Look up the index with nativeID
                            object idOrIndex = null;
                            if (sItr.Value.nativeID != null && sItr.Value.nativeID.Length > 0)
                            {
                                idOrIndex = sItr.Value.nativeID;
                                nativeIDList.Add(idOrIndex.ToString());
                            }
                        }//end if (boolcharge)
                    }//end foreach


            
            //get all the ids from pepxml
            //put them in dict, format =
            //ID,model,pep,mvh,mzf,tot,mat

            Dictionary<string, string> dict = new Dictionary<string, string>();
            try
            {
                using (XmlTextReader reader = new XmlTextReader(pepXMLFile))
                {
                    string peptideSeq = "";
                    string matchedIon = "";
                    string totalIon = "";
                    string mvh = "";
                    string mzFidelity = "";
                    string spectrumNativeID = "";
                    while (reader.Read())
                    {
                        if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("spectrum_query"))
                        {
                            spectrumNativeID = reader.GetAttribute("spectrumNativeID");
                            string assumed_charge = reader.GetAttribute("assumed_charge");

                            if (nativeIDList.Contains(spectrumNativeID) && assumed_charge == z.ToString())
                            {
                                XmlReader subReader = reader.ReadSubtree();
                                while (subReader.Read())
                                {
                                    if (subReader.NodeType.Equals(XmlNodeType.Element) && subReader.Name.Equals("search_hit"))
                                    {
                                        string hitRank = subReader.GetAttribute("hit_rank");
                                        if (hitRank == "1")
                                        {
                                            peptideSeq = subReader.GetAttribute("peptide");
                                            matchedIon = subReader.GetAttribute("num_matched_ions");
                                            totalIon = subReader.GetAttribute("tot_num_ions");
                                            XmlReader subReader2 = reader.ReadSubtree();
                                            while (subReader2.Read())
                                            {
                                                if (subReader2.NodeType.Equals(XmlNodeType.Element) && subReader.Name.Equals("search_score"))
                                                {
                                                    string name = subReader2.GetAttribute("name");
                                                    if (name == "mvh")
                                                    {
                                                        mvh = subReader2.GetAttribute("value");
                                                    }
                                                    if (name == "mzFidelity")
                                                    {
                                                        mzFidelity = subReader2.GetAttribute("value");
                                                    }
                                                }
                                            }//sub reader 2
                                        }
                                        break;
                                    }
                                }//sub reader 1

                                

                                //find the ions informations here. it is a long block
                                /////////////////////////////////////////////////////
                                //          long block
                                /////////////////////////////////////////////////////
                                string classStr = "";
                                
                                int len = peptideSeq.Length;
                                int spectrumIndex = sl.find(spectrumNativeID);
                                Spectrum spec1 = sl.spectrum(spectrumIndex, true);
                                MZIntensityPairList peaks = new MZIntensityPairList();
                                spec1.getMZIntensityPairs(ref peaks);
                                Set<Peak> peakList = new Set<Peak>();

                                //get TIC, and get the cutoff value, let TicCutoffPercentage=0.98
                                double TICValue = 0;
                                CVParamList list = spec1.cvParams;
                                foreach (CVParam CVP in list)
                                {
                                    if (CVP.name == "total ion current") TICValue = CVP.value;
                                }
                                double cutOffTIC = TicCutoffPercentage * TICValue;

                                Peptide peptide = new Peptide(peptideSeq, ModificationParsing.ModificationParsing_Auto, ModificationDelimiter.ModificationDelimiter_Brackets);
                                Fragmentation fragmentation = peptide.fragmentation(true, true);

                                //find the peaks, and find the peak cutoff; peak classes
                                //peak cutoff value is the value that: if intensity > cutoff, then upper class; if <= then lower class. 
                                double classABCutOff = 0;
                                double classBCCutOff = 0;
                                double intenThreshhold = 0;
                                int totalIntenClass = 0;

                                int indexIntenClassA = 0;
                                int indexIntenClassB = 0;
                                double[] intenArray = new double[peaks.Count];
                                //used during the foreach loop
                                int indexPeaks = 0;

                                //get all the peaks no matter how small they are
                                //then calculate the threshhold TIC
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

                                if (totalIntenClass % 7 == 0) indexIntenClassA = totalIntenClass / 7;
                                else indexIntenClassA = totalIntenClass / 7 + 1;
                                indexIntenClassB = indexIntenClassA * 2 + indexIntenClassA;
                                classABCutOff = intenArray[indexIntenClassA - 1];
                                classBCCutOff = intenArray[indexIntenClassB - 1];
                                //test here.
                                //Console.WriteLine("**********************");
                                Console.WriteLine("nativeID: " + spectrumNativeID);
                                //Console.WriteLine("test of intensity class:");
                                //Console.WriteLine("Tic: " + TICValue);
                                Console.WriteLine("class A/B cutoff: " + classABCutOff);
                                Console.WriteLine("class B/C cutoff: " + classBCCutOff);
                                Console.WriteLine("class C/D cutoff: " + intenThreshhold);
                                //Console.WriteLine();
                                //Console.WriteLine("intensity array: ");
                                //foreach (double dd in intenArray)
                                //{
                                //Console.Write(dd + "___");
                                //}


                                /////////////////////////////////////////////
                                //nm mode
                                /////////////////////////////////////////////
                                //count the total number of matched/unmatched ions in different classes
                                int nmionClassA = 0;
                                int nmionClassB = 0;
                                int nmionClassC = 0;
                                int nmionClassD = 0;
                                int nmion_matched = 0;
                                int nmion_tot_theretical = 0;

                                //define a list that contains fragment ions information
                                List<double> nmFragmentList = new List<double>();
                                for (int k = 1; k < len; k++)
                                {
                                    double[] bCharge = new double[z + 1];
                                    double[] yCharge = new double[z + 1];
                                    double[] bIntensity = new double[z + 1];
                                    double[] yIntensity = new double[z + 1];

                                    //return the b ion charge 1 Intensity, if matched
                                    for (int i = 1; i < z; i++)
                                    {
                                        //tags
                                        bCharge[i] = fragmentation.b(k, i);
                                        yCharge[i] = fragmentation.y(len - k, i);
                                        nmFragmentList.Add(bCharge[i]);
                                        nmFragmentList.Add(yCharge[i]);
                                    }
                                }


                                //the following will be addressing the match problem
                                double[] temp = new double[100];
                                int kkk = 0;
                                foreach (double mz in nmFragmentList)
                                {
                                    Peak matched = findNear(peakList, mz, 0.5);

                                    if (matched != null)
                                    {
                                        double intensity = matched.rankOrIntensity;
                                        if (intensity >= classABCutOff) nmionClassA++;
                                        else if (intensity >= classBCCutOff) nmionClassB++;
                                        else nmionClassC++;
                                        //test
                                        Console.WriteLine("mz matched is: " + matched.mz + "++++" + mz);
                                        temp[kkk] = matched.mz;
                                        kkk++;
                                    }
                                    else nmionClassD++;
                                   

                                    //test here
                                    //Console.WriteLine("intensity matched is: " + matched.rankOrIntensity);
                                    //Console.WriteLine("class info: " + nmionClassA + "," + nmionClassB + "," + nmionClassC + "," + nmionClassD);
                                }
                                Array.Sort(temp);
                                foreach (double dddd in temp) if (dddd > 0) Console.WriteLine(dddd);

                                nmion_matched = nmionClassA + nmionClassB + nmionClassC;
                                nmion_tot_theretical = nmFragmentList.Count;

                                string nmions = nmionClassA + "," + nmionClassB + "," + nmionClassC + "," + nmionClassD;

                                //end of nm model;


                                /////////////////////////////////////////////
                                //pm mode
                                /////////////////////////////////////////////
                                int pmionClassA = 0;
                                int pmionClassB = 0;
                                int pmionClassC = 0;
                                int pmionClassD = 0;
                                int pmion_matched = 0;
                                int pmion_tot_theretical = 0;

                                //define a list that contains fragment ions information
                                List<double> pmFragmentList = new List<double>();

                                int seqLength = peptideSeq.Length;
                                char[] seq = peptideSeq.ToCharArray();
                                int totalR = 0, totalK = 0, totalH = 0;

                                for (int i = 0; i < seqLength; ++i)
                                {
                                    if (seq[i] == 'R')
                                        ++totalR;
                                    else if (seq[i] == 'K')
                                        ++totalK;
                                    else if (seq[i] == 'H')
                                        ++totalH;
                                }

                                for (int c = 1; c < len; c++)
                                {
                                    int totalNR = 0, totalNK = 0, totalNH = 0, totalNL = 0;
                                    int totalCR = 0, totalCK = 0, totalCH = 0, totalCL = 0;
                                    for (int i = 0; i < c; ++i)
                                    {
                                        if (seq[i] == 'R') { ++totalNR; }
                                        else if (seq[i] == 'K') { ++totalNK; }
                                        else if (seq[i] == 'H') { ++totalNH; }
                                    }
                                    totalNL = c;
                                    totalCR = totalR - totalNR;
                                    totalCK = totalK - totalNK;
                                    totalCH = totalH - totalNH;
                                    totalCL = seqLength - totalNL;
                                    int R = totalNR - totalCR;
                                    int K = totalNK - totalCK;
                                    int H = totalNH - totalCH;
                                    int L = totalNL - totalCL;

                                    //tags
                                    string bion = peptideSeq.Substring(0, c);
                                    string yion = peptideSeq.Substring(c, len - c);

                                    double[] bCharge = new double[z + 1];
                                    double[] yCharge = new double[z + 1];



                                    if (z == 2)
                                    {
                                        double logit = totalNR * 2.7718555 + totalNK * 2.10656514 + totalNH * 2.21958398 + totalNL * 0.53519512 - totalCR * 1.37828528 - totalCK * 1.02145116 - totalCH * 1.0670024 - totalCL * 0.33021434;

                                        //for0/+2 situation, label = -2
                                        if (logit < -4.5087758)
                                        {
                                            yCharge[2] = fragmentation.y(len - c, 2);
                                            pmFragmentList.Add(yCharge[2]);

                                        }
                                        //for 0/+2, +1/+1
                                        else if (logit < -3.0816934)
                                        {
                                            bCharge[1] = fragmentation.b(c, 1);
                                            yCharge[1] = fragmentation.y(len - c, 1);
                                            yCharge[2] = fragmentation.y(len - c, 2);
                                            pmFragmentList.Add(bCharge[1]);
                                            pmFragmentList.Add(yCharge[1]);
                                            pmFragmentList.Add(yCharge[2]);
                                        }
                                        //for +1/+1
                                        else if (logit < 8.450989)
                                        {
                                            bCharge[1] = fragmentation.b(c, 1);
                                            yCharge[1] = fragmentation.y(len - c, 1);
                                            pmFragmentList.Add(bCharge[1]);
                                            pmFragmentList.Add(yCharge[1]);
                                        }
                                        //for +1/+1 +2/0
                                        else if (logit < 9.6695634)
                                        {
                                            bCharge[1] = fragmentation.b(c, 1);
                                            bCharge[2] = fragmentation.b(c, 2);
                                            yCharge[1] = fragmentation.y(len - c, 1);
                                            pmFragmentList.Add(bCharge[1]);
                                            pmFragmentList.Add(bCharge[2]);
                                            pmFragmentList.Add(yCharge[1]);
                                        }
                                        //for +2/0
                                        else
                                        {
                                            yCharge[2] = fragmentation.y(len - c, 2);
                                            pmFragmentList.Add(yCharge[2]);
                                        }
                                    } //end of charge 2

                                    if (z == 3)
                                    {
                                        double logit = totalNR * 1.99133766 + totalNK * 1.83742276 + totalNH * 1.88793026 + totalNL * 0.30088108 - totalCR * 1.15699882 - totalCK * 0.96586726 - totalCH * 1.07280076 - totalCL * 0.2461845;

                                        //for 0/+3 situation
                                        if (logit < -7.80125196)
                                        {
                                            yCharge[3] = fragmentation.y(len - c, 3);
                                            pmFragmentList.Add(yCharge[3]);

                                        }
                                        //for 0/+3, +1/+2
                                        else if (logit < -6.39099666)
                                        {
                                            bCharge[1] = fragmentation.b(c, 1);
                                            yCharge[2] = fragmentation.y(len - c, 2);
                                            yCharge[3] = fragmentation.y(len - c, 3);
                                            pmFragmentList.Add(bCharge[1]);
                                            pmFragmentList.Add(yCharge[2]);
                                            pmFragmentList.Add(yCharge[3]);
                                        }
                                        //for +1/+2
                                        else if (logit < 0.59062958)
                                        {
                                            bCharge[1] = fragmentation.b(c, 1);
                                            yCharge[2] = fragmentation.y(len - c, 2);
                                            pmFragmentList.Add(bCharge[1]);
                                            pmFragmentList.Add(yCharge[2]);
                                        }
                                        //for +1/+2 +2/+1
                                        else if (logit < 2.0965649)
                                        {
                                            bCharge[1] = fragmentation.b(c, 1);
                                            bCharge[2] = fragmentation.b(c, 2);
                                            yCharge[1] = fragmentation.y(len - c, 1);
                                            yCharge[2] = fragmentation.y(len - c, 2);
                                            pmFragmentList.Add(bCharge[1]);
                                            pmFragmentList.Add(bCharge[2]);
                                            pmFragmentList.Add(yCharge[1]);
                                            pmFragmentList.Add(yCharge[2]);
                                        }
                                        //for +2/+1
                                        else if (logit < 9.7471158)
                                        {
                                            bCharge[2] = fragmentation.b(c, 2);
                                            yCharge[1] = fragmentation.y(len - c, 1);
                                            pmFragmentList.Add(bCharge[2]);
                                            pmFragmentList.Add(yCharge[1]);
                                        }
                                        //for +2/+1 +3/0
                                        else if (logit < 10.9444849)
                                        {
                                            bCharge[2] = fragmentation.b(c, 2);
                                            bCharge[3] = fragmentation.b(c, 3);
                                            yCharge[1] = fragmentation.y(len - c, 1);
                                            pmFragmentList.Add(bCharge[2]);
                                            pmFragmentList.Add(bCharge[3]);
                                            pmFragmentList.Add(yCharge[1]);
                                        }
                                        //for +3/0
                                        else
                                        {
                                            bCharge[3] = fragmentation.b(c, 3);
                                            pmFragmentList.Add(bCharge[3]);
                                        }
                                    } //end of charge 3
                                }//end of peptide

                                //the following will be addressing the match problem
                                foreach (double mz in pmFragmentList)
                                {
                                    Peak matched = findNear(peakList, mz, 0.5);

                                    if (matched != null)
                                    {
                                        double intensity = matched.rankOrIntensity;
                                        if (intensity >= classABCutOff) pmionClassA++;
                                        else if (intensity >= classBCCutOff) pmionClassB++;
                                        else pmionClassC++;
                                    }
                                    else pmionClassD++;
                                }

                                pmion_matched = pmionClassA + pmionClassB + pmionClassC;
                                pmion_tot_theretical = pmFragmentList.Count;

                                string pmions = pmionClassA + "," + pmionClassB + "," + pmionClassC + "," + pmionClassD;
                                classStr = nmions + "," + pmions;
                                ///////////////////////////////////////////////
                                //  end of the long block
                                //////////////////////////////////////////////

                                string ionString = spectrumNativeID + "," + model + "," + peptideSeq + "," + mvh + "," + mzFidelity + "," + totalIon + "," + matchedIon + "," + classStr;
                                dict.Add(spectrumNativeID, ionString);
                                Console.WriteLine("spectrumNativeID added...");
                                Console.WriteLine(dict[spectrumNativeID]);
                            }
                        }
                    }
                }
            }
            catch (Exception exc)
            {
                //in order to avoid exception, I need to put the test.idpXML in \obj\debug folder. wonder why?
                Console.Write("\r\nError in reading pepXML file, please check IDPicker configuration and try again\r\n");
                Console.Write("Close");
                throw new Exception(exc.Message);
            }

            return dict;
        }

        /// ///////////////////////////////////////////////////////////////////
        /// ///////////////////////////////////////////////////////////////////
        /// end of task
        /// ///////////////////////////////////////////////////////////////////
        /// ///////////////////////////////////////////////////////////////////
        


        public static void Main(string[] args)
        {
            /*
            //for idpReaderPlus
            string idpXMLFile_nm = "H:\\home\\dwang\\complex12-10-10\\test\\test_nm.idpXML";
            string idpXMLFile_pm = "H:\\home\\dwang\\complex12-10-10\\test\\test_pm.idpXML";
            string pepXMLFile_nm = "H:\\home\\dwang\\complex12-10-10\\test\\test_nm.pepXML";
            string pepXMLFile_pm = "H:\\home\\dwang\\complex12-10-10\\test\\test_pm.pepXML";

            string output = "H:\\home\\dwang\\complex12-10-10\\test\\output.csv";

            string mzMLFile = "H:\\home\\dwang\\complex12-10-10\\test\\test.mzML";
            string model_nm = "nm";
            string model_pm = "pm";
            double ltq = 0.95;
            double orbi = 0.98;
            int z = 3;

            Dictionary<string, string> dict = pepXMLReader(idpXMLFile_nm, pepXMLFile_nm, mzMLFile, model_nm, 0.95, 3);

            

            Console.WriteLine("writing into csv file...");
            using (StreamWriter file = new StreamWriter(output))
            {
                string head = "ID,model,pep,mvh,mzF,mat,tot,A,B,C,D,A,B,C,D";
                file.WriteLine(head);
                foreach (string nativeID in dict.Keys)
                {
                    file.WriteLine(dict[nativeID]);
                }
            }
             */ 
             

            //original program
            //idpReader   
            
            if (args.Length != 2 )
            {
                Console.WriteLine("\nError!");
                Console.WriteLine("\nUsage: idpReader.exe <idpXML file> <corresponding mzML file> <cutoffpercentage> <precursor charge>");
                return;
            }
            if (!File.Exists(args[0]))
            {
                Console.Write("\r\nError: Cannot find input idpXML file!\n");
                return;
            }
            if (!File.Exists(args[1]))
            {
                Console.Write("\r\nError: Cannot find corresponding mzML file!\n");
                return;
            }
           
            //string idpXMLFile_nm = "H:\\home\\dwang\\complex12-10-10\\test\\test_nm.idpXML";
            //string mzMLFile = "H:\\home\\dwang\\complex12-10-10\\test\\test.mzML";
            idpReader(args[0], args[1], 0.95, 3);
             
        }

    }//end IdpReader
}
