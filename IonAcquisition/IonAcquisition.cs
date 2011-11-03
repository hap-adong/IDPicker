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


namespace IonRetrieval
{
    class IonAcquisition
    {
        //grab idpxml confident ids
        //analyze it by MM nm model
        //analyze it by MM Basophile model
        //make comparison
        public static Dictionary<string, string> ionAcqusition(string idpXMLFile, string pepXMLFile, string mzMLFile, double TicCutoffPercentage, int z)
        {
            List<string> nmmzIntensityList = new List<string>();
            List<string> pmmzIntensityList = new List<string>();
            List<string> mzIntensityList = new List<string>();
            nmmzIntensityList.Add("id,mz,intensity");
            pmmzIntensityList.Add("id,mz,intensity");
            mzIntensityList.Add("id,mz,intensity");

            //list number of: matched ion,unique ion, replciate ion
            List<string> testList = new List<string>();

            //first get IDs from idpXML file
            Dictionary<string, string> idpxmlNativeIDDict = new Dictionary<string, string>();
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
                        bool boolCharge = sItr.Value.id.charge.Equals(z);
                        if (boolCharge)
                        {
                            string rawPepSequence = vi.ToString();
                            string interpretation = vi.ToSimpleString();
<<<<<<< HEAD
=======
                            Console.WriteLine(interpretation);
                            Console.WriteLine(vi.peptide.sequence);
>>>>>>> 80fc9e47b4dafd28bcc96f478ae26543036b9ff6
                            // Look up the index with nativeID
                            object idOrIndex = null;
                            if (sItr.Value.nativeID != null && sItr.Value.nativeID.Length > 0)
                            {
                                idOrIndex = sItr.Value.nativeID;
                                idpxmlNativeIDDict.Add(idOrIndex.ToString(), interpretation);
                            }
                        }//end if (boolcharge)
                    }//end foreach



            //get all the ids from pepxml
            //put them in dict, format =
            //ID,model,pep,mvh,mzf,tot,mat

            Dictionary<string, string> dict = new Dictionary<string, string>();

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

                        if (idpxmlNativeIDDict.ContainsKey(spectrumNativeID) && assumed_charge == z.ToString())
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

                            Peptide peptide = new Peptide(idpxmlNativeIDDict[spectrumNativeID], ModificationParsing.ModificationParsing_Auto, ModificationDelimiter.ModificationDelimiter_Brackets);
                            Fragmentation fragmentation = peptide.fragmentation(true, true);

                            //find the peaks, and find the peak cutoff; peak classes
                            //peak cutoff value is the value that: if intensity > cutoff, then upper class; if <= then lower class. 
                            double classABCutOff = 0;
                            double classBCCutOff = 0;
                            double intenThreshhold = 0;

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

                            //new method of judging the threshold
                            //code from MM.
                            double relativeIntensity = 0.0f;
                            for (int i = 0; relativeIntensity < TicCutoffPercentage && i < peaks.Count; i++)
                            {

                                relativeIntensity += (intenArray[i] / TICValue);
                                intenThreshhold = intenArray[i];
                            }

                            //old method.
                            //if currTIC>=cutoff, then break
                            /*
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
                            */

                            //then based on that, generate a new peaklist that contains only ABC peaks
                            //then calculate the intensity classes
                            foreach (MZIntensityPair mzIntensity in peaks)
                            {
                                if (mzIntensity.intensity >= intenThreshhold)
                                {
                                    mzIntensityList.Add(spectrumNativeID + "," + mzIntensity.mz + "," + mzIntensity.intensity);
                                    Peak p = new Peak(mzIntensity.mz, mzIntensity.intensity);
                                    peakList.Add(p);
                                }
                            }
                            double peakDataSize = peakList.Count;
                            indexIntenClassA = (int)(Math.Round(peakDataSize / 7));
                            indexIntenClassB = indexIntenClassA + (int)(Math.Round(peakDataSize * 2 / 7));
                            classABCutOff = intenArray[indexIntenClassA - 1];
                            classBCCutOff = intenArray[indexIntenClassB - 1];

                            /////////////////////////////////////////////
                            //nm mode
                            /////////////////////////////////////////////
                            
                            //count the total number of matched/unmatched ions in different classes
                            int nmionClassA = 0;
                            int nmionClassB = 0;
                            int nmionClassC = 0;
                            int nmionClassD = 0;

                            //define a list that contains fragment ions information
                            List<double> nmFragmentList = new List<double>();
                            //defines a list that contain the matched.mz and then find out the replicates
                            List<double> matchedmzList = new List<double>();
                            for (int k = 1; k < len; k++)
                            {
                                double[] bCharge = new double[z + 1];
                                double[] yCharge = new double[z + 1];
                                double[] bIntensity = new double[z + 1];
                                double[] yIntensity = new double[z + 1];

                                //return the b ion charge 1 Intensity, if matched
                                //important here!!!!!!
                                //for MM nm model, no +3 is modeled. 
                                for (int i = 1; i < z; i++)
                                {
                                    //tags
                                    bCharge[i] = fragmentation.b(k, i);
                                    yCharge[i] = fragmentation.y(len - k, i);
                                    nmFragmentList.Add(bCharge[i]);
                                    nmFragmentList.Add(yCharge[i]);
                                }
                            }


                            
                            //could delete the sort here.
                            nmFragmentList.Sort();
                            
                            foreach (double mz in nmFragmentList)
                            {
                                Peak matched = Package.findClose(peakList, mz, 0.5);
                                if (matched != null)
                                {
                                    double intensity = matched.rankOrIntensity;
                                    nmmzIntensityList.Add(spectrumNativeID + "," + matched.mz + "," + intensity);
                                    //matchedmzList.Add(matched.mz);
                                    if (intensity >= classABCutOff) nmionClassA++;
                                    else if (intensity >= classBCCutOff) nmionClassB++;
                                    else nmionClassC++;
                                }
                                else
                                {
                                    nmionClassD++;
                                    nmmzIntensityList.Add(spectrumNativeID + "," + mz + ",0");
                                }
                            }

                            //find out how many ions are identically matched. 
                            //Console.WriteLine(spectrumNativeID + "...........");
                            //int numBefore = matchedmzList.Count;
                            //List<double> replicate = Package.findCommon(matchedmzList);

                            //testList.Add(numBefore + "," + replicate.Count);
                           
                            //nmion_matched = nmionClassA + nmionClassB + nmionClassC;
                            //nmion_tot_theretical = nmFragmentList.Count;
                            //string nmions = nmionClassA + "," + nmionClassB + "," + nmionClassC + "," + nmionClassD;
                            //end of nm model;
                            

                            
                                
                            /////////////////////////////////////////////
                            //pm mode
                            /////////////////////////////////////////////
                            
                            int pmionClassA = 0;
                            int pmionClassB = 0;
                            int pmionClassC = 0;
                            int pmionClassD = 0;

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
                                    double logit = totalNR*0.5249139 + totalNK*0.3574616 + totalNH*0.3863936 + totalNL*0.2388976 - totalCR*0.4124967 - totalCK*0.3049158 -totalCH*0.2443213 - totalCL*0.2337753;

                                    //for 0/+3 situation
                                    if (logit < -6.386273)
                                    {
                                        yCharge[3] = fragmentation.y(len - c, 3);
                                        pmFragmentList.Add(yCharge[3]);

                                    }
                                    //for 0/+3, +1/+2
                                    else if (logit < -3.287443)
                                    {
                                        bCharge[1] = fragmentation.b(c, 1);
                                        yCharge[2] = fragmentation.y(len - c, 2);
                                        yCharge[3] = fragmentation.y(len - c, 3);
                                        pmFragmentList.Add(bCharge[1]);
                                        pmFragmentList.Add(yCharge[2]);
                                        pmFragmentList.Add(yCharge[3]);
                                    }
                                    //for +1/+2
                                    else if (logit < -1.222116)
                                    {
                                        bCharge[1] = fragmentation.b(c, 1);
                                        yCharge[2] = fragmentation.y(len - c, 2);
                                        pmFragmentList.Add(bCharge[1]);
                                        pmFragmentList.Add(yCharge[2]);
                                    }
                                    //for +1/+2 +2/+1
                                    else if (logit < 1.390020)
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
                                    else if (logit < 3.624574)
                                    {
                                        bCharge[2] = fragmentation.b(c, 2);
                                        yCharge[1] = fragmentation.y(len - c, 1);
                                        pmFragmentList.Add(bCharge[2]);
                                        pmFragmentList.Add(yCharge[1]);
                                    }
                                    //for +2/+1 +3/0
                                    else if (logit < 7.117430)
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
                                Peak matched = Package.findClose(peakList, mz, 0.5);
                                if (matched != null)
                                {
                                    double intensity = matched.rankOrIntensity;
                                    pmmzIntensityList.Add(spectrumNativeID + "," + matched.mz + "," + intensity);
                                    matchedmzList.Add(matched.mz);
                                    if (intensity >= classABCutOff) pmionClassA++;
                                    else if (intensity >= classBCCutOff) pmionClassB++;
                                    else pmionClassC++;
                                }
                                else
                                { 
                                    pmionClassD++;
                                    pmmzIntensityList.Add(spectrumNativeID + "," + mz + ",0");
                                }
                            }

                            //find out how many ions are identically matched. 
                            int numBefore = matchedmzList.Count;
                            List<double> replicate = Package.findCommon(matchedmzList);

                            testList.Add(numBefore + "," + replicate.Count);

                            ///////////////////////////////////////////////
                            //  end of the long block
<<<<<<< HEAD
                            //////////////////////////////////////////////
=======
                            ///////////////////////////////////////////////
>>>>>>> 80fc9e47b4dafd28bcc96f478ae26543036b9ff6

                            //string ionString = spectrumNativeID + "," + peptideSeq + "," + mvh + "," + mzFidelity + "," + totalIon + "," + matchedIon + "," + nmions;
                            //dict.Add(spectrumNativeID, ionString);
                        }//end of if (containKey...)



                    }//end of if(Elemetary.Element and spectrum querry)
                }//end of whileRead
            }//end of using file
            StreamWriter file_nmmzIntensity = new StreamWriter("C:\\Documents and Settings\\wangd5\\Model12-28-10\\test\\nmmzIntensity.csv");
            StreamWriter file_pmmzIntensity = new StreamWriter("C:\\Documents and Settings\\wangd5\\Model12-28-10\\test\\pmmzIntensity.csv");
            StreamWriter file_mzIntensity = new StreamWriter("C:\\Documents and Settings\\wangd5\\Model12-28-10\\test\\mzIntensity.csv");
            StreamWriter file_test = new StreamWriter("C:\\Documents and Settings\\wangd5\\Model12-28-10\\test\\test.csv");
            foreach (string ss in nmmzIntensityList) file_nmmzIntensity.WriteLine(ss);
            foreach (string ss in pmmzIntensityList) file_pmmzIntensity.WriteLine(ss);
            foreach (string ss in mzIntensityList) file_mzIntensity.WriteLine(ss);
            foreach (string ss in testList) file_test.WriteLine(ss);
            return dict;
        }

        static void Main(string[] args)
        {
<<<<<<< HEAD
            //string idpxml_nm = "C:\\Documents and Settings\\wangd5\\Model12-28-10\\yeast_ORBI_nm.idpXML";
            //string idpxml_pm = "C:\\Documents and Settings\\wangd5\\Model12-28-10\\idpXML\\yeast_ORBI_pm.idpXML";
            //string mzML = "C:\\Documents and Settings\\wangd5\\Model12-28-10\\yeast_ORBI.mzML";
            //IonRetrieval.NativeIDComparison.nativeIDComp(idpxml_nm, idpxml_pm, mzML);
            //Console.WriteLine();

            //the following for pepxmlreader...........
            Console.WriteLine("usage: ionAcqusition <idpxml> <pepxml> <mzml>");
            string fileName = Path.GetFileNameWithoutExtension(args[0]);
            string filePath = Path.GetDirectoryName(args[0]);
            string csvFile = Path.Combine(filePath, fileName) + "_" + 3.ToString() + ".csv";

            Dictionary<string, string> dict = ionAcqusition(args[0], args[1], args[2], 0.98, 3);
            Console.WriteLine("writing into csv file...");
            using (StreamWriter file = new StreamWriter(csvFile))
            {
                string head = "ID,pepSeq,mvh,mzF,totIon,matchedIon,A,B,C,D";
                file.WriteLine(head);
                foreach (string nativeID in dict.Keys)
                {
                    file.WriteLine(dict[nativeID]);
                }
            }
=======
            //string idpxml_nm = "Z:\\home\\dwang\\fragmentation\\LTQ\\yeast20090403\\MVH\\binary+2model\\mam_20090403x_Yeast_60nguL_1_naive.idpXML";
            //string idpxml_pm = "Z:\\home\\dwang\\fragmentation\\LTQ\\yeast20090403\\MVH\\basophilenew\\mam_20090403x_Yeast_60nguL_1.idpXML";
            //string mzML = "Z:\\home\\dwang\\fragmentation\\LTQ\\yeast20090403\\mam_20090403x_Yeast_60nguL_1.mzXML";
            string idpXMLFile = "Z:\\home\\dwang\\fragmentation\\LTQ\\mam_20090403x_Yeast_60nguL_1.idpXML";
            string pepXMLFile = "Z:\\home\\dwang\\fragmentation\\LTQ\\mam_20090403x_Yeast_60nguL_1.pepXML";
            int z=2;
            string output = "Z:\\home\\dwang\\fragmentation\\LTQ\\yeast20090403\\MVH\\charge2baso.csv";

            //IDPicker.Workspace workspace = new IDPicker.Workspace();
            //Package.loadWorkspace(ref workspace, idpXMLFile);


            //MSDataFile foo = new MSDataFile(mzMLFile);
            //SpectrumList sl = foo.run.spectrumList;

            //foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
            //    foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
            //        foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
            //        {
            //            IDPicker.ResultInstance ri = sItr.Value.results[1];
            //            IDPicker.VariantInfo vi = ri.info.peptides.Min;
            //            bool boolCharge = sItr.Value.id.charge.Equals(z);
            //            if (boolCharge)
            //            {
            //                string rawPepSequence = vi.ToString();
            //                string interpretation = vi.ToSimpleString();
            //                Console.WriteLine(interpretation);
            //                Console.WriteLine(vi.peptide.sequence);
            //                // Look up the index with nativeID
            //                object idOrIndex = null;
            //                if (sItr.Value.nativeID != null && sItr.Value.nativeID.Length > 0)
            //                {
            //                    idOrIndex = sItr.Value.nativeID;
            //                    idpxmlNativeIDDict.Add(idOrIndex.ToString(), interpretation);
            //                }
            //            }//end if (boolcharge)
            //        }//end foreach


            List<string> list = new List<string>();
            using (XmlTextReader reader = new XmlTextReader(pepXMLFile))
            {
                string matchedIon = "";
                string totalIon = "";
                string spectrumNativeID = "";
                while (reader.Read())
                {
                    if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("spectrum_query"))
                    {
                        spectrumNativeID = reader.GetAttribute("spectrumNativeID");
                        string assumed_charge = reader.GetAttribute("assumed_charge");

                        if (/*idpxmlNativeIDDict.ContainsKey(spectrumNativeID) && */assumed_charge == z.ToString())
                        {
                            XmlReader subReader = reader.ReadSubtree();
                            while (subReader.Read())
                            {
                                if (subReader.NodeType.Equals(XmlNodeType.Element) && subReader.Name.Equals("search_hit"))
                                {
                                    string hitRank = subReader.GetAttribute("hit_rank");
                                    if (hitRank == "1")
                                    {
                                        matchedIon = subReader.GetAttribute("num_matched_ions");
                                        totalIon = subReader.GetAttribute("tot_num_ions");
                                        list.Add(matchedIon + "," + totalIon);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            
            TextWriter file = new StreamWriter(output);
            foreach (string ss in list)
                file.WriteLine(ss);
            file.Close();

                 
            //IonRetrieval.NativeIDComparison.nativeIDComp(idpxml_nm, idpxml_pm, mzML);
            //Console.WriteLine();

            //the following for pepxmlreader...........
            //Console.WriteLine("usage: ionAcqusition <idpxml> <pepxml> <mzml>");
            //string fileName = Path.GetFileNameWithoutExtension(args[0]);
            //string filePath = Path.GetDirectoryName(args[0]);
            //string csvFile = Path.Combine(filePath, fileName) + "_" + 3.ToString() + ".csv";

            //Dictionary<string, string> dict = ionAcqusition(args[0], args[1], args[2], 0.98, 3);
            //Console.WriteLine("writing into csv file...");
            //using (StreamWriter file = new StreamWriter(csvFile))
            //{
            //    string head = "ID,pepSeq,mvh,mzF,totIon,matchedIon,A,B,C,D";
            //    file.WriteLine(head);
            //    foreach (string nativeID in dict.Keys)
            //    {
            //        file.WriteLine(dict[nativeID]);
            //    }
            //}
>>>>>>> 80fc9e47b4dafd28bcc96f478ae26543036b9ff6
        }
    }
}
