using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml;
using System.IO;
using System.Data;
using pwiz.CLI.proteome;



namespace Test
{
    class Program
    {
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



        static void Main(string[] args)
        {
            float x = 79.9f;
            if (Math.Round(x) == 80)
                Console.WriteLine("kick");



            //some job
            List<string> listA = new List<string>();
            List<string> listB = new List<string>();

            string groupA = "C:\\Users\\Dong\\Desktop\\groupA.txt";
            TextReader txA = new StreamReader(groupA);
            string line;
            while (( line = txA.ReadLine()) != null)
            {
                //Console.WriteLine(line);
                listA.Add(line);  
            }
            txA.Close();


            string groupB = "C:\\Users\\Dong\\Desktop\\groupB.txt";
            TextReader txB = new StreamReader(groupB);
            while ((line = txB.ReadLine()) != null)
            {
                System.Text.RegularExpressions.Regex m = new System.Text.RegularExpressions.Regex(", ");

                string[] strs = m.Split(line);
                for (int i = 0; i < strs.Length; i++)
                {
                    string str = strs[i];
                    //Console.WriteLine(str);
                    listB.Add(str);
                }
            }
            txB.Close();

            Console.WriteLine("group A: " + listA.Count);
            Console.WriteLine("group B: " + listB.Count);
            //to make comparison
            List<string> A = new List<string>();
            List<string> B = new List<string>();
            A = Package.findCommon(listA);
            B = Package.findCommon(listB);
            Console.WriteLine("group A normalized: " + A.Count);
            Console.WriteLine("group B normalized: " + B.Count);

            int common =0;
            foreach (string geneA in A)
            {
                foreach (string geneB in B)
                {
                    if (geneB == geneA)
                    {
                        common++;
                    }
                }
            }
            Console.WriteLine("number of common genes: " + common);

            //read the database
            //List<string> proteinList = readDatabase("X:\\wangd5\\s1_20_seq-reverse.fasta");

            //get all confident peptides
            //List<string> nativeIDList = new List<string>();
            //string xml = "K:\\project vandy\\dataset\\chymo\\naive\\mam_20090622n_KChen_Yeast_Chymotrypsin_184min_1.idpXML";
            //IDPicker.Workspace workspace = new IDPicker.Workspace();
            //Package.loadWorkspace(ref workspace, xml);
            //foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
            //    foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
            //        foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
            //        {
            //            IDPicker.ResultInstance ri = sItr.Value.results[1];
            //            IDPicker.VariantInfo vi = ri.info.peptides.Min;
            //            string nativeID = sItr.Value.nativeID;
            //            nativeIDList.Add(nativeID);
            //            Console.WriteLine(nativeID);
            //        }


            //Console.WriteLine("reading pepxml file...");
            //List<string> list = new List<string>();
            //string pepXMLFile = "X:\\wangd5\\pepxml\\myrimatch\\mam_121007n_RKO_200ug_IEF_A01.pepXML";
            //using (XmlTextReader reader = new XmlTextReader(pepXMLFile))
            //{
            //    while (reader.Read())
            //    {
            //        if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("spectrum_query"))
            //        {
            //            string spectrumNativeID = reader.GetAttribute("spectrumNativeID");
            //            string combo = "";

            //            if (nativeIDList.Contains(spectrumNativeID))
            //            {
            //                XmlReader subReader = reader.ReadSubtree();
            //                while (subReader.Read())
            //                {
            //                    if (subReader.NodeType.Equals(XmlNodeType.Element) && subReader.Name.Equals("search_hit"))
            //                    {
            //                        //string hitRank = subReader.GetAttribute("hit_rank");
                                    
            //                        string peptideSeq = "";
            //                        string rna = "0";
            //                        string mvh = "";
            //                        string mzfidelity = "";

            //                        peptideSeq = subReader.GetAttribute("peptide");
                                    
            //                        foreach (var protein in proteinList)
            //                        {
            //                            if (protein.Contains(peptideSeq))
            //                            {
            //                                rna = "1";
            //                                break;
            //                            }
            //                        }

            //                        XmlReader subReader2 = reader.ReadSubtree();
            //                        while (subReader2.Read())
            //                        {
            //                            if (subReader2.NodeType.Equals(XmlNodeType.Element) && subReader.Name.Equals("search_score"))
            //                            {
            //                                string name = subReader2.GetAttribute("name");
            //                                if (name == "mvh")
            //                                {
            //                                    mvh = subReader2.GetAttribute("value");
            //                                }
            //                                if (name == "mzFidelity")
            //                                {
            //                                    mzfidelity = subReader2.GetAttribute("value");
            //                                }
            //                            }
            //                        }

            //                        combo += "," + mvh + "," + mzfidelity + "," + rna;

            //                    }//end of node
            //                }//sub reader 1
            //                if (combo != "")
            //                    list.Add(combo);
            //            }//end of if
                        
            //        }
            //    }
            //}//end of using

            //Console.WriteLine("writing csv file...");
            //string output = "X:\\wangd5\\idpXML_FDR1.00\\score evaluation\\confranks.csv";
            //TextWriter tw = new StreamWriter(output);
            //foreach (string ss in list)
            //    tw.WriteLine(ss);

        }//end of main

    }
}
