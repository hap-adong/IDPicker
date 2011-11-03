using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using pwiz.CLI.proteome;
using System.Text.RegularExpressions;


namespace UPS
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

        //find replicate in a double list
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

        //find replicate in a string list
        public static List<string> findCommon(List<string> inputList)
        {
            Dictionary<string, int> uniqueStore = new Dictionary<string, int>();
            List<string> finalList = new List<string>();

            foreach (string currValue in inputList)
            {
                if (!uniqueStore.ContainsKey(currValue))
                {
                    uniqueStore.Add(currValue, 0);
                }
                else finalList.Add(currValue);
            }
            return finalList;
        }
        
        //find common ids shared bytwo string list
        //input should be no replicates, or else must perform removeDuplicate
        public static List<string> findCommon(List<string> inputList1, List<string> inputList2)
        {
            List<string> finalList = new List<string>();

            foreach (string currValue in inputList1)
            {
                if (inputList2.Contains(currValue))
                {
                    finalList.Add(currValue);
                }
            }
            return finalList;
        }

        //find unique string list from a string list
        public static List<string> removeDuplicate(List<string> inputList)
        {
            Dictionary<string, int> uniqueStore = new Dictionary<string, int>();
            List<string> finalList = new List<string>();

            foreach (string currValue in inputList)
            {
                if (!uniqueStore.ContainsKey(currValue))
                {
                    uniqueStore.Add(currValue, 0);
                    finalList.Add(currValue);
                }
            }
            return finalList;
        }

        //read a fasta database
        //remove the reverse sequences
        //write forward sequences into a proteinList
        public static List<string> readDatabase(string database)
        {
            List<string> proteinList = new List<string>();
            ProteomeDataFile foo = new ProteomeDataFile(database);
            ProteinList p1 = foo.proteinList;

            for (int i = 0; i < p1.size(); i++)
            {
                Protein pro = p1.protein(i, true);
                if (!pro.id.Contains("rev"))
                {
                    string sequence = p1.protein(i, true).sequence;
                    proteinList.Add(sequence);
                }

            }
            return proteinList;
        }

        //given a assembly.xml file, grab peptide ids based on charge state.
        //z=0, output all the peptides
        public static List<string> PepSecurity(string xml, int z)
        {
            List<string> peptideList = new List<string>();
            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, xml);
            foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
                foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
                    foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
                    {
                        IDPicker.ResultInstance ri = sItr.Value.results[1];
                        IDPicker.VariantInfo vi = ri.info.peptides.Min;
                        string rawPepSequence = vi.ToString();
                        string pepSequence = vi.peptide.sequence;
                        if (z == 0)
                        {
                            if (sItr.Value.id.charge.Equals(2) || sItr.Value.id.charge.Equals(3) || sItr.Value.id.charge.Equals(4))
                            {
                                peptideList.Add(pepSequence);
                            }
                        }
                       
                        //if (z == 0) peptideList.Add(pepSequence);
                        else
                        {
                            if (sItr.Value.id.charge.Equals(z)) peptideList.Add(pepSequence);
                        }
                    }
            List<string> uniPepList = Package.removeDuplicate(peptideList);

            return uniPepList;
        }


        //search tree
        //input a peptide list that do not contain replicates
        //make a single peptide tree and map to proteinList
        public static int pepMaping(List<string> peptideList, List<string> proteinList)
        {
            int match = 0;
            foreach (var peptide in peptideList)
            {
                //given a single peptide, make a single tree.
                foreach (var protein in proteinList)
                {
                    if (protein.Contains(peptide))
                    {
                        match++;
                        break;
                    }
                }
            }
            return match;
        }


        //read xml files to get peptide, non-redundant
        //read fasta to get non-reversed protein
        //make them into list, respectively
        //note the charge state in pepSecurity
        public static void pep2RNASeq()
        {
            string MM = "X:\\wangd5\\idpXML_FDR1.00\\myrimatch\\Assemble_MM.xml";
            string XT = "X:\\wangd5\\idpXML_FDR1.00\\tandem\\Assemble_XT.xml";
            string SQ = "X:\\wangd5\\idpXML_FDR1.00\\sequest\\Assemble_SQ.xml";

            //string MM = "X:\\wangd5\\idpXML_FDR0.05\\myrimatch\\Assemble_MM.xml";
            //string XT = "X:\\wangd5\\idpXML_FDR0.05\\tandem\\output_all\\Assemble_XT.xml";
            //string SQ = "X:\\wangd5\\idpXML_FDR0.05\\sequest\\Assemble_SQ.xml";
            List<string> MM_pep = Package.PepSecurity(MM, 0);
            List<string> XT_pep = Package.PepSecurity(XT, 0);
            List<string> SQ_pep = Package.PepSecurity(SQ, 0);
           
            //try ways to get overlap information.
            List<string> MM_XT_SQ_pep = new List<string>();
            
            foreach (string pep in MM_pep) MM_XT_SQ_pep.Add(pep);
            foreach (string pep in XT_pep) MM_XT_SQ_pep.Add(pep);
            foreach (string pep in SQ_pep) MM_XT_SQ_pep.Add(pep);
            MM_XT_SQ_pep = Package.removeDuplicate(MM_XT_SQ_pep);

            //define 7 different regions of venn diagram
            //MM_ex, XT_ex, SQ_ex
            //MM_XT, MM_SQ, XT_SQ
            //MM_XT_SQ
            List<string> MM_ex = new List<string>();
            List<string> XT_ex = new List<string>();
            List<string> SQ_ex = new List<string>();

            List<string> MM_XT = new List<string>();
            List<string> MM_SQ = new List<string>();
            List<string> XT_SQ = new List<string>();

            List<string> MM_XT_SQ = new List<string>();

            foreach (string pep in MM_XT_SQ_pep)
            {
                if (MM_pep.Contains(pep) && !XT_pep.Contains(pep) && !SQ_pep.Contains(pep)) MM_ex.Add(pep);
                if (XT_pep.Contains(pep) && !MM_pep.Contains(pep) && !SQ_pep.Contains(pep)) XT_ex.Add(pep);
                if (SQ_pep.Contains(pep) && !XT_pep.Contains(pep) && !MM_pep.Contains(pep)) SQ_ex.Add(pep);

                if (MM_pep.Contains(pep) && XT_pep.Contains(pep) && !SQ_pep.Contains(pep)) MM_XT.Add(pep);
                if (MM_pep.Contains(pep) && !XT_pep.Contains(pep) && SQ_pep.Contains(pep)) MM_SQ.Add(pep);
                if (!MM_pep.Contains(pep) && XT_pep.Contains(pep) && SQ_pep.Contains(pep)) XT_SQ.Add(pep);

                if (MM_pep.Contains(pep) && XT_pep.Contains(pep) && SQ_pep.Contains(pep)) MM_XT_SQ.Add(pep);
            }

            //output overlap information.
            Console.WriteLine("==========output overlap information==========");
            Console.WriteLine("%%7 different regions of venn diagram");
            Console.WriteLine("%%MM_ex, XT_ex, SQ_ex, MM_XT, MM_SQ, XT_SQ, MM_XT_SQ");
            Console.WriteLine(MM_ex.Count + "," + XT_ex.Count + "," + SQ_ex.Count + "||" + MM_XT.Count + "," + MM_SQ.Count + "," + XT_SQ.Count + "||" + MM_XT_SQ.Count);
            Console.WriteLine();


            string fasta = "X:\\wangd5\\s1_microarray95_seq-reverse.fasta";
            List<string> pro = Package.readDatabase(fasta);
            int mm_ex = Package.pepMaping(MM_ex, pro);
            int xt_ex = Package.pepMaping(XT_ex, pro);
            int sq_ex = Package.pepMaping(SQ_ex, pro);
            int mm_xt = Package.pepMaping(MM_XT, pro);
            int mm_sq = Package.pepMaping(MM_SQ, pro);
            int xt_sq = Package.pepMaping(XT_SQ, pro);
            int mm_xt_sq = Package.pepMaping(MM_XT_SQ, pro);
            Console.WriteLine("==========output mapping information==========");
            Console.WriteLine("%%MM_ex, XT_ex, SQ_ex, MM_XT, MM_SQ, XT_SQ, MM_XT_SQ");
            Console.WriteLine(mm_ex + "," + xt_ex + "," + sq_ex + "||" + mm_xt + "," + mm_sq + "," + xt_sq + "||" + mm_xt_sq);
        }

        //the speccharge2score get each spectrum-peptide, and then map it to scores. 
        public static Dictionary<string, double> SpecCharge2Score(string idpxml, string qonversion, List<string> list)
        {
            ///<summary>
            ///given qonversion and idpxml files, find the scores, DecoyState, making a dictionary. 
            ///make it peptide-eccentric. 
            ///</summary>
            
            //get peptide sequence
            string name = Path.GetFileNameWithoutExtension(idpxml);
            

            Dictionary<string, string> dic_spectrum = new Dictionary<string, string>(); 

            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, idpxml);
            int repeatedIndex = 0;
            foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
                foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
                    foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
                    {
                        IDPicker.ResultInstance ri = sItr.Value.results[1];
                        IDPicker.VariantInfo vi = ri.info.peptides.Min;
                        string rawPepSequence = vi.ToString();
                        string pepSequence = vi.peptide.sequence;
                        string index = sItr.Value.id.index.ToString();
                        int z = sItr.Value.id.charge;
                        string id = name + "." + index + "," + z.ToString();
                        if (!dic_spectrum.ContainsKey(id))
                        {
                            dic_spectrum.Add(id, pepSequence);
                        }
                        else repeatedIndex++;
                    }

            //Console.WriteLine("repeated index (comes from spectrum with ambiguous charge state) is: " + repeatedIndex);

            //then read qonverion.txt
            //get index, TotalScore, DecoyState
            //add pep and RNASeq_support information

            //dic_peptide is what we wanna return
            Dictionary<string, double> finaldic = new Dictionary<string, double>();
            string t;
            int unmatchedPep = 0;
            TextReader file = new StreamReader(qonversion);
            while ((t = file.ReadLine()) != null)
            {
                if (t.Contains("("))
                {
                    Regex r = new Regex(" +");
                    string[] str = r.Split(t);
                    string Index = str[2];
                    double TotalScore = Convert.ToDouble(str[8]);
                    string charge = str[3];
                    string DecoyState = str[4];
                    string id = name + "." + Index + "," + charge;
                    if (dic_spectrum.ContainsKey(id))
                    {
                        string pep = dic_spectrum[id];
                        string key = id + "," + pep + "," + DecoyState;
                        finaldic.Add(key, TotalScore);
                        list.Add(key);
                    }
                    else unmatchedPep++;

                }
            }
            //Console.WriteLine("unmatched peptide (peptide that was not found in the csv file) is: " + unmatchedPep);
            Console.WriteLine("--one set xml-qonversion done with unmatched (unexpected): " + unmatchedPep);
            return finaldic;

        }//end of spectra2RNASeq


        //the Pep2Score get each peptide, and then map it to highest scores. 
        public static Dictionary<string, double> Pep2Score(string idpxml, string qonversion, List<string> list)
        {
            ///<summary>
            ///given qonversion and idpxml files, find the the highest scores, DecoyState, charge, N/C specificity
            ///make it peptide-eccentric. 
            ///</summary>

            //get peptide sequence
            string name = Path.GetFileNameWithoutExtension(idpxml);


            Dictionary<string, string> dic_pep = new Dictionary<string, string>();
            Dictionary<string, string> dic_spec = new Dictionary<string, string>();
            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, idpxml);
            foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
                foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
                    foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
                    {
                        IDPicker.ResultInstance ri = sItr.Value.results[1];
                        IDPicker.VariantInfo vi = ri.info.peptides.Min;
                        string rawPepSequence = vi.ToString();
                        string pepSequence = vi.peptide.sequence;
                        string index = sItr.Value.id.index.ToString();
                        int z = sItr.Value.id.charge;
                        string id = pepSequence + "," + z.ToString();
                        dic_spec.Add(name + "," + index + "," + z, pepSequence);
                    }

            //Console.WriteLine("repeated index (comes from spectrum with ambiguous charge state) is: " + repeatedIndex);

            //then read qonverion.txt
            //get index, TotalScore, DecoyState
            //add pep and RNASeq_support information

            //dic_peptide is what we wanna return
            Dictionary<string, double> final = new Dictionary<string, double>();
            string t;
            TextReader file = new StreamReader(qonversion);
            while ((t = file.ReadLine()) != null)
            {
                if (t.Contains("("))
                {
                    Regex r = new Regex(" +");
                    string[] str = r.Split(t);
                    string Index = str[2];
                    double TotalScore = Convert.ToDouble(str[8]);
                    string charge = str[3];
                    string DecoyState = str[4];
                    string id = name + "," + Index + "," + charge;
                    if (dic_spec.ContainsKey(id))
                    {
                        string pep = dic_spec[id];
                        string key = pep + "," + DecoyState;
                        list.Add(key);
                        if (final.ContainsKey(key))
                        {
                            if (TotalScore > final[key])
                            {
                                final[key] = TotalScore;
                            }
                        }
                        else
                        {
                            final.Add(key, TotalScore);
                        }
                        
                    }

                }
            }
            //Console.WriteLine("unmatched peptide (peptide that was not found in the csv file) is: " + unmatchedPep);
            Console.WriteLine("one set of xml and qonversion files done");
            return final;

        }//end of spectra2RNASeq


        //specific to RKO dataset   or SW480 dataset
        //spectra based analysis
        public static Dictionary<string, double> UPS(string path, List<string> list)
        { 
            Dictionary<string, double> dic_spectrum2score = new Dictionary<string, double>();
            //keep the keys of the dictionary into a list

            //keep track of those repetition, in which spectrumID + charge + peptide repeats (likely to be 0)
            int repeats = 0;
            
            string xml = path  + ".idpXML";
            string qonversion = path + "-qonversion.txt";
            Dictionary<string, double> dic = Package.SpecCharge2Score(xml, qonversion, list);
            foreach (string key in dic.Keys)
            {
                if (!dic_spectrum2score.ContainsKey(key))
                {
                    dic_spectrum2score.Add(key, dic[key]);
                }
                else
                {
                    repeats++;
                }
            }
            Console.WriteLine("the final repeats in this search engine is: " + repeats + "\n");
            return dic_spectrum2score;
        }

    }//end of package.
}//end of namespace
