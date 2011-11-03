using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data;
using System.IO;
namespace UPS
{
    class Program
    {
        static void Main(string[] args)
        {
            List<string> ProteinList = Package.readDatabase("Z:\\home\\dwang\\fragmentation\\UPS\\Sigma49-reverse.fasta");

            string path_n = "Z:\\home\\dwang\\fragmentation\\UPS\\naive\\klc_031308p_cptac_study6_6_QC1.idpXML";
            string path_b = "Z:\\home\\dwang\\fragmentation\\UPS\\MVH\\basophilenew\\klc_031308p_cptac_study6_6_QC1.idpXML";

            for (int z = 0; z <= 4; z++)
            {
                int num_n = 0, num_b = 0, num_c = 0;
                Console.WriteLine("starting analyzing z===" + z);
                List<string> pep_n = Package.PepSecurity(path_n, z);
                List<string> pep_b = Package.PepSecurity(path_b, z);
                List<string> pep_common = Package.findCommon(pep_n, pep_b);
                foreach (string ss in pep_n)
                {
                    foreach (var protein in ProteinList)
                    {
                        if (protein.Contains(ss))
                        {
                            num_n++;
                            break;
                        }
                    }
                }
                foreach (string ss in pep_b)
                {
                    foreach (var protein in ProteinList)
                    {
                        if (protein.Contains(ss))
                        {
                            num_b++;
                            break;
                        }
                    }
                }
                foreach (string ss in pep_common)
                {
                    foreach (var protein in ProteinList)
                    {
                        if (protein.Contains(ss))
                        {
                            num_c++;
                            break;
                        }
                    }
                }

                Console.WriteLine("number of naive: " + num_n + "/" + pep_n.Count);
                Console.WriteLine("number of baso: " + num_b + "/" + pep_b.Count);
                Console.WriteLine("number of common: " + num_c + "/" + pep_common.Count);
                Console.WriteLine();
            }
            
            
            /*

            //===========================================================
            string test_mm = "X:\\wangd5\\UPS\\MM\\UPS_FDR0.05.idpXML";
            string test_sq = "X:\\wangd5\\UPS\\SQ\\UPS_FDR0.05.idpXML";
            string test_xt = "X:\\wangd5\\UPS\\XT\\UPS_FDR0.05.idpXML";
            string test_p = "X:\\wangd5\\UPS\\evaluation\\p.csv";

            //get peptides in p
            List<string> p = new List<string>();
            TextReader file_p = new StreamReader(test_p);
            DataTable dt = CSV.CsvParser.Parse(file_p, true);
            foreach (DataRow dr in dt.Rows)
            {
                string pep = dr[2].ToString();
                p.Add(pep);
            }

            List<string> m = Package.PepSecurity(test_mm, 0);
            List<string> x = Package.PepSecurity(test_xt, 0);
            List<string> s = Package.PepSecurity(test_sq, 0);

            m = Package.removeDuplicate(m);
            x = Package.removeDuplicate(x);
            s = Package.removeDuplicate(s);
            p = Package.removeDuplicate(p);

            List<string> pm = Package.findCommon(p, m);



            List<string> mx = Package.findCommon(m, x);
            List<string> ms = Package.findCommon(m, s);
            List<string> xs = Package.findCommon(x, s);
            List<string> mxs = Package.findCommon(mx, s);




            ///////////////////////////////////////////////////////////////////////////////////////

            //read the fdr1.0.csv, to get all the peptide sequences that identified by 3 engines\
            //no matter how low the score is
            //also, the file contains the RNASeq information
            Dictionary<string, string> dic_RNASeq = new Dictionary<string, string>();
            TextReader file_csv = new StreamReader("X:\\wangd5\\UPS\\evaluation\\PepRNA.csv");
            //TextReader file_csv = new StreamReader("X:\\wangd5\\idpXML_FDR1.00\\score evaluation\\fdr1.0.csvn");
            DataTable table = CSV.CsvParser.Parse(file_csv, true);
            foreach (DataRow dr in table.Rows)
            {
                string pep = dr[0].ToString();
                string RNASeq = dr[1].ToString();
                dic_RNASeq.Add(pep, RNASeq);
            }


            int m_n = 0;
            int x_n = 0;
            int s_n = 0;
            int p_n = 0;
            int ms_n = 0;
            int mx_n = 0;
            int xs_n = 0;
            int mxs_n = 0;
            int pm_n = 0;

            foreach (string ss in m)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") m_n++;
                }
                else Console.WriteLine("crap");

            }
            foreach (string ss in x)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") x_n++;
                }

            }

            foreach (string ss in s)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") s_n++;
                }

            }

            foreach (string ss in p)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") p_n++;
                }

            }

            foreach (string ss in ms)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") ms_n++;
                }

            }
            foreach (string ss in mx)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") mx_n++;
                }

            }
            foreach (string ss in xs)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") xs_n++;
                }

            }
            foreach (string ss in mxs)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") mxs_n++;
                }

            }

            foreach (string ss in pm)
            {
                if (dic_RNASeq.ContainsKey(ss))
                {
                    if (dic_RNASeq[ss] == "1") pm_n++;
                }

            }

            Console.WriteLine("m: " + m.Count + "=" + m_n);
            Console.WriteLine("x: " + x.Count + "=" + x_n);
            Console.WriteLine("s: " + s.Count + "=" + s_n);
            Console.WriteLine("p: " + p.Count + "=" + p_n);

            Console.WriteLine("mx: " + mx.Count + "=" + mx_n);
            Console.WriteLine("ms: " + ms.Count + "=" + ms_n);
            Console.WriteLine("xs: " + xs.Count + "=" + xs_n);

            Console.WriteLine("mxs: " + mxs.Count + "=" + mxs_n);

            Console.WriteLine("pm: " + pm.Count + "=" + pm_n);
            
             */
        }
    }
}
