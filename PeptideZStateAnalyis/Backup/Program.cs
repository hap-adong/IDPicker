using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml;
using System.ComponentModel;
using System.IO;
using System.Data;
using pwiz.CLI;
using pwiz.CLI.data;
using pwiz.CLI.msdata;
using pwiz.CLI.proteome;

namespace ConsoleApplication1
{
    
    /// <summary>
    /// A peak in a mass spectrum.
    /// </summary>
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


    class FragZStateAnalyzer
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

        //define a int var to record the number of confident peptides (different from csv report)
        public static int numberSpectrumIdentified;
        //dictionary contain ion index and m/z
        public static Dictionary<int, double> ionDict = new Dictionary<int, double>();

        private static T getAttributeAs<T>(XmlTextReader reader, string attribute, bool throwIfAbsent)
        {
            if (reader.MoveToAttribute(attribute))
            {
                TypeConverter c = TypeDescriptor.GetConverter(typeof(T));
                if (c == null || !c.CanConvertFrom(typeof(string)))
                    throw new Exception("unable to convert from string to " + typeof(T).Name);
                T value = (T)c.ConvertFromString(reader.Value);
                reader.MoveToElement();
                return value;
            }
            else if (throwIfAbsent)
                throw new Exception("missing required attribute \"" + attribute + "\"");
            else if (typeof(T) == typeof(string))
                return (T)TypeDescriptor.GetConverter(typeof(T)).ConvertFromString(String.Empty);
            else
                return default(T);
        }

        // read idpxml, extract spectra id.charge, save to a dictionary
        public static Dictionary<int, string> getAllIdCharge(string path)
        {
            //record two dic, idtScanDict record the "number of Id, ID+charge", and finalScanDict record the "number of id, pepNum+ID+charge"
            Dictionary<int, string> idtScanDict = new Dictionary<int, string>();

            //File.Exists might be very special. the default folder should be sth like \bin\debug
            if (!File.Exists(path))
            {
                Console.Write("\r\nError: Cannot find idpXML file!");
                Console.Write("\r\nPlease make sure to input the full path and name of the file\n");
                return idtScanDict;
            }
            Console.WriteLine("\nparsing to get the charge +3 dictionary...\n");
            try
            {
                using (XmlTextReader reader = new XmlTextReader(path))
                {
                    while (reader.Read())
                    {
                        if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("spectrum"))
                        {
                            // Read the spectrum tag
                            //  <spectrum id="614" nativeID="614" index="196" z="1" mass="569.32" 
                            //   time="16.7" targets="82" decoys="0" results="1">
                            numberSpectrumIdentified++;

                            string nativeID = getAttributeAs<string>(reader, "id", false);
                            Match m = Regex.Match(nativeID, @"scan=(\d+)");

                            int z = getAttributeAs<int>(reader, "z", true);
                            //int index = getAttributeAs<int>(reader, "index", true);
                            string idtScan = nativeID + "." + Convert.ToString(z);
                            idtScanDict.Add(numberSpectrumIdentified, idtScan);

                        }

                        //here, there is an exception thrown out by some point of the dictionary
                        if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("id"))
                        {
                            int pepSequenceNum = getAttributeAs<int>(reader, "peptide", false);
                            string idtScan = Convert.ToString(pepSequenceNum) + "," + idtScanDict[numberSpectrumIdentified];
                            idtScanDict[numberSpectrumIdentified] = idtScan;

                        }
                    }
                }

                return idtScanDict;
            }
            catch (Exception exc)
            {
                //in order to avoid exception, I need to put the test.idpXML in \obj\debug folder. wonder why?
                Console.Write("\r\nError in reading idpXML file, please check IDPicker configuration and try again\r\n");
                Console.Write("Close");
                throw new Exception(exc.Message);
            }
        }

        public static Dictionary<int, string> replaceSequenceWithId(Dictionary<int, string> dict, string path)
        {
            //what idpxml <peptide> like:
            // <peptide id="218" sequence="ERFGFSGFPVTEDGK" mass="1671.7892855281002" unique="1" NTerminusIsSpecific="1" CTerminusIsSpecific="1">
            // <locus id="370" offset="141" />
            //</peptide>


            //get pepID from dict (first part of it, limited by ",")
            List<int> list = new List<int>(dict.Keys);
            foreach (int n in list)
            {
                int peptideID;
                Match m = (Regex.Match(dict[n], @"(\d+),"));
                if (m.Success)
                {
                    peptideID = Convert.ToInt32(m.Groups[1].Value);
                    string[] arr = dict[n].Split(',');
                    int lengthArr = arr.Length;


                    //read the idpxml file
                    try
                    {
                        using (XmlTextReader reader = new XmlTextReader(path))
                        {
                            while (reader.Read())
                            {
                                if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("peptide"))
                                {
                                    int pepID = getAttributeAs<int>(reader, "id", false);
                                    if (peptideID == pepID)
                                    {
                                        string pepSequence = getAttributeAs<string>(reader, "sequence", false);
                                        dict[n] = pepSequence + "," + arr[lengthArr - 1];

                                    }

                                }
                            }
                        }

                    }
                    catch (Exception exc)
                    {
                        Console.Write("\r\nError in getting sequence (method getSequence)\r\n");
                        Console.Write("Close");
                        throw new Exception(exc.Message);
                    }


                }
            }


            return dict;
        }


        //get all the +3 information
        public static Dictionary<int, string> getCharge3(Dictionary<int, string> dict)
        {
            Dictionary<int, string> charge3Dict = new Dictionary<int, string>();
            List<int> list = new List<int>(dict.Keys);
            int i = 0;
            foreach (int n in list)
            {
                int charge;

                Match m = (Regex.Match(dict[n], @"\.(\d+)"));
                if (m.Success)
                {

                    charge = Convert.ToInt32(m.Groups[1].Value);
                    if (charge == 3)
                    {
                        i++;
                        charge3Dict[i] = dict[n];
                    }
                }
            }
            return charge3Dict;
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
        public static string parsePepSequence(string pepSequence)
        {
            int numR = 0;
            int numK = 0;
            int numH = 0;
            int numAA = 0;


            char[] charArr = pepSequence.ToCharArray();
            int len = pepSequence.Length;
            string AA = "";


            foreach (char c in charArr)
            {
                if (c == 'R') numR++;
                if (c == 'K') numK++;
                if (c == 'H') numH++;
            }
            numAA = len - numR - numK - numH;
            AA = numR.ToString() + "," + numK.ToString() + "," + numH.ToString() + "," + numAA.ToString();
            return AA;


        }

            
        public static List<string> parseTable(Dictionary<int, string> dict, string path)
        {
            List<string> seqIntensityList = new List<string>();
            seqIntensityList.Add("nativeID,peptideSequence,NR,NK,NH,NAA,CR,CK,CH,CAA,b1Intensity,b2Intensity,b3Intensity,y1Intensity,y2Intensity,y3Intensity,basePeak,TIC");

            MSDataFile foo = new MSDataFile(path);
            SpectrumList sl = foo.run.spectrumList;

            if (!File.Exists(path))
            {
                Console.Write("\r\nError: Cannot find mzML file!");
                Console.Write("\r\nPlease make sure to input the full path and name of the file\n");
                Console.Write("parsing stopped!");
                return seqIntensityList;
            }

            Console.WriteLine("\nparsing to get the final table...\n");
            foreach (int n in dict.Keys)
            {
                Console.WriteLine("peptide ->...");
                //get the predicted peptide sequence and nativeID
                string[] s = dict[n].Split(',');
                string[] ss = s[1].Split('.');
                string pepSequence = s[0];
                string nativeID = ss[0];
                int len = pepSequence.Length;

                
                //get parameter 1: spectrum information->peakList

                int spectrumIndex = sl.find(nativeID);
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
                    if (CVP.name == "base peak m/z")
                    {
                        basePeakValue = CVP.value;
                    }
                    if (CVP.name == "total ion current")
                    {
                        TICValue = CVP.value;
                    }
                }
                double cutOffBasePeakValue = 0.05 * basePeakValue;
                string basePeak = basePeakValue.ToString();
                string TIC = TICValue.ToString();

                Peptide peptide = new Peptide(pepSequence);
                Fragmentation fragmentation = peptide.fragmentation(true, true);

                Console.WriteLine("prepare peaklist ->...");
                foreach (MZIntensityPair pp in peaks)
                {
                    Peak p = new Peak(pp.mz, pp.intensity);
                    peakList.Add(p);
                }

                //call the method
                Console.WriteLine("fragment ion ->...");
                for (int k = 1; k < len; k++)
                {
                    double b1Intensity = 0;
                    double b2Intensity = 0;
                    double b3Intensity = 0;
                    double y1Intensity = 0;
                    double y2Intensity = 0;
                    double y3Intensity = 0;

                    string bion = pepSequence.Substring(0, k);
                    string yion = pepSequence.Substring(k, len - k);
                    string bionAANumber = parsePepSequence(bion);
                    string yionAANumber = parsePepSequence(yion);


                    double bCharge1 = fragmentation.b(k, 1);
                    double bCharge2 = fragmentation.b(k, 2);
                    double bCharge3 = fragmentation.b(k, 3);
                    double yCharge1 = fragmentation.y(len - k, 1);
                    double yCharge2 = fragmentation.y(len - k, 2);
                    double yCharge3 = fragmentation.y(len - k, 3);

                    //return the b ion charge 1 Intensity, if matched
                    Peak matchedb1 = findNear(peakList, bCharge1, 0.5);
                    if (matchedb1 != null)
                    {
                        b1Intensity = matchedb1.rankOrIntensity;
                        if (b1Intensity < cutOffBasePeakValue)
                        {
                            b1Intensity = 0;
                        }
                    }

                    //return the b ion charge 2 Intensity, if matched
                    Peak matchedb2 = findNear(peakList, bCharge2, 0.5);
                    if (matchedb2 != null)
                    {
                        b2Intensity = matchedb2.rankOrIntensity;
                        if (b2Intensity < cutOffBasePeakValue)
                        {
                            b2Intensity = 0;
                        }
                    }

                    //return the b ion charge 3 Intensity, if matched
                    Peak matchedb3 = findNear(peakList, bCharge3, 0.5);
                    if (matchedb3 != null)
                    {
                        b3Intensity = matchedb3.rankOrIntensity;
                        if (b3Intensity < cutOffBasePeakValue)
                        {
                            b3Intensity = 0;
                        }
                    }

                    //return the y ion charge 1 Intensity, if matched
                    Peak matchedy1 = findNear(peakList, yCharge1, 0.5);
                    if (matchedy1 != null)
                    {
                        y1Intensity = matchedy1.rankOrIntensity;
                        if (y1Intensity < cutOffBasePeakValue)
                        {
                            y1Intensity = 0;
                        }
                    }

                    //return the y ion charge 2 Intensity, if matched
                    Peak matchedy2 = findNear(peakList, yCharge2, 0.5);
                    if (matchedy2 != null)
                    {
                        y2Intensity = matchedy2.rankOrIntensity;
                        if (y2Intensity < cutOffBasePeakValue)
                        {
                            y2Intensity = 0;
                        }
                    }

                    //return the y ion charge 3 Intensity, if matched
                    Peak matchedy3 = findNear(peakList, yCharge3, 0.5);
                    if (matchedy3 != null)
                    {
                        y3Intensity = matchedy3.rankOrIntensity;
                        if (y3Intensity < cutOffBasePeakValue)
                        {
                            y3Intensity = 0;
                        }
                    }


                    string ionString = b1Intensity.ToString() + "," + b2Intensity.ToString() + "," + b3Intensity.ToString() + "," + y1Intensity.ToString() + "," + y2Intensity.ToString() + "," + y3Intensity.ToString();
                    if (ionString != "0,0,0,0,0,0")
                    {
                        string finalString = nativeID + "," + pepSequence + "," + bionAANumber + "," + yionAANumber + "," + ionString + "," + basePeak + "," + TIC;
                        seqIntensityList.Add(finalString);
                        Console.WriteLine("seqIntensityList successfully added ->...");
                    }

                }
                Console.WriteLine("\n\n");

            }
            //parse done.
            Console.WriteLine("\nparse getTable done\n");
            return seqIntensityList;
        }



        static void Main(string[] args)
        {
            if (args.Length < 3)
            {
                Console.WriteLine("\nError!");
                Console.WriteLine("\nUsage: getTable.exe <idpXML file> <corresponding mzML file> <output csv file>");
                return;
            }
            if (!File.Exists(args[2]))
            {
                Console.Write("\r\nError: Cannot find output csv file!\n");
                return;
            }

            FragZStateAnalyzer analyzer = new FragZStateAnalyzer();
            analyzer.readWorkspace(args[0]);

            
            Dictionary<int, string> getDict = getAllIdCharge(args[0]);
            Dictionary<int, string> get2Dict = replaceSequenceWithId(getDict, args[0]);
            Dictionary<int, string> tempDict = getCharge3(get2Dict);

            List<string> finalList = parseTable(tempDict, args[1]);

            //write tableList into a csv file. 
            using (StreamWriter file = new StreamWriter(args[2]))
            {
                Console.WriteLine("\nStarting writing into csv file...\n");
                foreach (string line in finalList) file.WriteLine(line);
            }

            Console.WriteLine("\nwriting complete\n");
        }
    }
}
