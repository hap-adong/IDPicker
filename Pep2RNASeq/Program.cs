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
using System.Data;
using pwiz.CLI.analysis;
using pwiz.CLI.data;
using pwiz.CLI.msdata;
using pwiz.CLI.proteome;
using SourceList = System.Collections.Generic.Set<IDPicker.SourceInfo>;

namespace Pep2RNASeq
{
    class temp
    {
        static void Main(string[] args)
        {
            
            //===========================================================
            string test_mm = "X:\\wangd5\\SW480\\MM\\FDR0.05\\mm.xml";
            string test_sq = "X:\\wangd5\\SW480\\SQ\\FDR0.05\\sq.xml";
            string test_xt = "X:\\wangd5\\SW480\\XT\\FDR0.05\\xt.xml";
            string test_p = "X:\\wangd5\\SW480\\evaluation\\p.csv";

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
            TextReader file_csv = new StreamReader("X:\\wangd5\\SW480\\evaluation\\PepRNA_FDR1.csv");
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
            


            string path_mm = "X:\\wangd5\\SW480\\MM\\FDR1.00\\mam_012808n_SW480_200ug_";
            string path_sq = "X:\\wangd5\\SW480\\SQ\\FDR1.00\\mam_012808n_SW480_200ug_";
            string path_xt = "X:\\wangd5\\SW480\\XT\\FDR1.00\\mam_012808n_SW480_200ug_";
            
            List<string> keyslist = new List<string>();
            Dictionary<string, double> dic_mm = Package.RKO(path_mm, keyslist);
            Dictionary<string, double> dic_sq = Package.RKO(path_sq, keyslist);
            Dictionary<string, double> dic_xt = Package.RKO(path_xt, keyslist);

            keyslist = Package.removeDuplicate(keyslist);

            //merge the three dictionaries into one. 
            Dictionary<string, string> dic_merge = new Dictionary<string,string>();
            int misses = 0;
            foreach (string key in keyslist)
            {
                string mm = "";
                string sq = "";
                string xt = "";
                string[] str = key.Split(',');
                string pep = str[2];
                if (dic_RNASeq.ContainsKey(pep))
                {
                    string RNA = dic_RNASeq[pep];

                    if (dic_mm.ContainsKey(key))
                    {
                        mm = dic_mm[key].ToString();
                    }
                    else mm = "0";
                    if (dic_xt.ContainsKey(key))
                    {
                        xt = dic_xt[key].ToString();
                    }
                    else xt = "0";
                    if (dic_sq.ContainsKey(key))
                    {
                        sq = dic_sq[key].ToString();
                    }
                    else sq = "0";
                    dic_merge.Add(key, mm + "," + xt + "," + sq + "," + RNA);
                }
                else misses++;
                
            }

            Console.WriteLine("the number of pep-rna misses is: " + misses);
            //write into file
            string output = "X:\\wangd5\\SW480\\evaluation\\merge.csv";
            TextWriter file = new StreamWriter(output);
            file.WriteLine("spectrum,charge,pep,decoystate,mm,xt,sq,rna");
            foreach (string key in dic_merge.Keys)
            {
                file.WriteLine(key + "," + dic_merge[key]);
            }
            
             
        }
        
    }
}
