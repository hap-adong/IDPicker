using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Xml;

namespace pepXMLReader
{
    class Package
    {
        public static double pepXMLReader(string pepxml, string output)
        {
            double totalions_pep = 0;
            int totalions = 0;
            int totalpep = 0;
            Console.WriteLine("reading pepxml file...");
            List<string> list = new List<string>();
            using (XmlTextReader reader = new XmlTextReader(pepxml))
            {
                while (reader.Read())
                {
                    if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("spectrum_query"))
                    {
                        string spectrumNativeID = reader.GetAttribute("spectrumNativeID");
                        string charge = reader.GetAttribute("assumed_charge");
                        if (charge == "3")
                        {
                            //test here
                            //Console.WriteLine("spectrumID= " + spectrumNativeID);
                            XmlReader subReader = reader.ReadSubtree();
                            while (subReader.Read())
                            {
                                if (subReader.NodeType.Equals(XmlNodeType.Element) && subReader.Name.Equals("search_hit"))
                                {
                                    string hitRank = subReader.GetAttribute("hit_rank");
                                    string peptideSeq = subReader.GetAttribute("peptide");
                                    string num_matched_ions = subReader.GetAttribute("num_matched_ions");
                                    string tot_num_ions = subReader.GetAttribute("tot_num_ions");
                                    int total_num_ions = Convert.ToInt16(tot_num_ions);
                                    totalions += total_num_ions;
                                    //test
                                    //Console.WriteLine("hit_rank= " + hitRank);
                                    string mvh = "";
                                    string mzfidelity = "";
                                    string hgt = "";
                                    string rst = "";
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
                                                mzfidelity = subReader2.GetAttribute("value");
                                            }
                                            if (name == "hgt")
                                            {
                                                hgt = subReader2.GetAttribute("value");
                                            }
                                            if (name == "rst")
                                            {
                                                rst = subReader2.GetAttribute("value");
                                                //test
                                                //Console.WriteLine("rst= " + rst);
                                            }
                                        }
                                    }//end of subReader2
                                    totalpep++;
                                    string finalstring = pepxml + "," + num_matched_ions + "," + tot_num_ions + "," + mvh + "," + mzfidelity + "," + hgt + "," + rst;
                                    list.Add(finalstring);
                                }//end of node: search_hit
                            }//sub reader 1
                        }
                        
                    }//end of if(spectruem_querry)
                } //while read pepXML
            }//end of using

            Console.WriteLine("writing csv file...");
            TextWriter tw = new StreamWriter(output);
            tw.WriteLine("pepxml,machion,totions,mvh,mf,hgt,rst");
            foreach (string ss in list)
            {
                tw.WriteLine(ss);
            }
            tw.Flush();
            tw.Close();

            totalions_pep = totalions*1.0/totalpep;
            return totalions_pep;
        }//end of method: pepXMLReader
    }
}
