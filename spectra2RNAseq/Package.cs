using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using pwiz.CLI.proteome;
using System.Data;
using System.Text.RegularExpressions;


namespace Spectra2RNAseq
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

        
        public static Dictionary<string, string> spectra2RNASeq(string idpxml, string qonversion, Dictionary<string, string> dic_RNASeq)
        {
            ///<summary>
            ///now the problem is:
            ///d a idpQonvert file, read in all the peptides from pepxml without any fdr cutoff
            ///ead an idpxml file, read in all peptides that passed the fdr cutoff
            ///ow I restrict my rows to fdr-cutoff peptides, which might generate a limited pool
            ///ext is to consider fdr 0.1 cutoff, it will generate enough "0"s in the RNASeq_support column
            ///ifferent cutoff values here might have different merits. we'll find out which is the best
            ///</summary>
            ///next to check 1. whether 0.1 cutoff has the same idpqonvert file or not; Yes, the same idpqonvert file
            ///2. generate new idpxml with fdr cutoff=0.1 done. 
            ///file path: X:\wangd5\idpXML_FDR1.00\myrimatch\FDR0.1

            
            //then read idpxml file
            //get peptide sequence
            string name = Path.GetFileNameWithoutExtension(idpxml);

            Dictionary<string, string> dic_peptide = new Dictionary<string, string>();
            //dic_spectra is what we want to return
            Dictionary<string, string> dic_spectra = new Dictionary<string, string>();
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
                        string id = name + "." + index + "." + z.ToString() ;
                        if (!dic_peptide.ContainsKey(id))
                        {
                            bool Nspecificity = vi.peptide.NTerminusIsSpecific;
                            bool Cspecitificy = vi.peptide.CTerminusIsSpecific;
                            
                            string value = pepSequence + "," + Nspecificity.ToString() + "," + Cspecitificy.ToString();
                            dic_peptide.Add(id, pepSequence);
                            dic_spectra.Add(id, value);
                        }
                        else repeatedIndex++;
                    }

            Console.WriteLine("repeated index (comes from spectrum with ambiguous charge state) is: " + repeatedIndex);

            //then read qonverion.txt
            //get index, charges , TotalScore, DecoyState
            //add pep and RNASeq_support information
            string t;
            int unmatchedPep = 0;
            TextReader file = new StreamReader(qonversion);
            while ((t = file.ReadLine()) != null)
            {
                if (t.Contains("("))
                {
                    Regex r = new Regex(" +");
                    string[] str = r.Split(t);
                    string NativeID = str[1];
                    string Index = str[2];
                    string charge = str[3];
                    string DecoyState = str[4];
                    string TotalScore = str[8];
                    string FDR = str[9];
                    string id = name + "." + Index + "." + charge;
                    if (dic_peptide.ContainsKey(id))
                    {
                        string pep = dic_peptide[id];
                        if (dic_RNASeq.ContainsKey(pep))
                        {
                            string RNASeq = dic_RNASeq[pep];
                            dic_spectra[id] = NativeID + "," + dic_spectra[id] + "," + charge + "," + DecoyState + "," + TotalScore + "," + FDR + "," + RNASeq;
                        }
                        else unmatchedPep++;
                        
                    }
                }
            }
            Console.WriteLine("unmatched peptide (peptide that was not found in the csv file) is: " + unmatchedPep);
            return dic_spectra;

        }//end of spectra2RNASeq


    }//end of package.
}//end of namespace
