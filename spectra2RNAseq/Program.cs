using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Data;
using System.Text.RegularExpressions;
using pwiz.CLI.msdata;


namespace Spectra2RNAseq
{
    class Program
    {
        static void Main(string[] args)
        {



            //TextReader ttt = new StreamReader("X:\\wangd5\\idpXML_FDR1.00\\myrimatch\\FDR0.1\\results\\pep.csv");
            //DataTable t1 = CSV.CsvParser.Parse(ttt, true);
            //List<string> temp = new List<string>();
            //foreach (DataRow dr in t1.Rows)
            //{
            //    string pep = dr[2].ToString();
            //    temp.Add(pep);
            //}
            //temp = Package.removeDuplicate(temp);
            //Console.WriteLine(temp.Count);




            //read previously completed file: X:\wangd5\idpXML_FDR1.00\myrimatch\FDR0.1\results\fdr1.0.csv
            //record each peptide sequence and outcome. 
            Dictionary<string, string> dic_RNASeq = new Dictionary<string, string>();
            TextReader file_csv = new StreamReader("X:\\wangd5\\idpXML_FDR1.00\\myrimatch\\FDR0.5\\results\\fdr1.0.csv");
            DataTable table = CSV.CsvParser.Parse(file_csv, true);
            foreach (DataRow dr in table.Rows)
            {
                string pep = dr[0].ToString();
                string outcome = dr[4].ToString();
                dic_RNASeq.Add(pep, outcome);
            }


            Dictionary<string, string> dic_all = new Dictionary<string, string>();
            Dictionary<string, string> dic_tmp = new Dictionary<string, string>();

            string path = "X:\\wangd5\\idpXML_FDR1.00\\myrimatch\\FDR0.5\\mam_121007n_RKO_200ug_IEF_";
            for (int i = 1; i <= 3; i++)
            {
                for (int j = 1; j <= 10; j++)
                {
                    string let = "";
                    string num = "";
                    switch (i)
                    {
                        case 1:
                            let = "A";
                            break;
                        case 2:
                            let = "B";
                            break;
                        case 3:
                            let = "C";
                            break;
                        default: 
                            Console.WriteLine("error");
                            break;
                    }
                    switch (j)
                    {
                        case 1:
                            num = "01";
                            break;
                        case 2:
                            num = "02";
                            break;
                        case 3:
                            num = "03";
                            break;
                        case 4:
                            num = "04";
                            break;
                        case 5:
                            num = "05";
                            break;
                        case 6:
                            num = "06";
                            break;
                        case 7:
                            num = "07";
                            break;
                        case 8:
                            num = "08";
                            break;
                        case 9:
                            num = "09";
                            break;
                        case 10:
                            num = "10";
                            break;
                        default:
                            Console.WriteLine("error");
                            break;
                    }
                    string ext = let + num;
                    string xml = path + ext + ".idpXML";
                    string qonversion = path + ext + "-qonversion.txt";
                    dic_tmp = Package.spectra2RNASeq(xml, qonversion, dic_RNASeq);
                    foreach (string key in dic_tmp.Keys)
                    {
                        dic_all.Add(key, dic_tmp[key]);
                    }
                }
            }
            

            string output = "X:\\wangd5\\idpXML_FDR1.00\\myrimatch\\FDR0.5\\results\\raw.csv";
            TextWriter file = new StreamWriter(output);
            file.WriteLine("id,nativeID,pep,Nspec,Cspec,charge,DecoyState,TotalScore,FDR,RNASeq");
            foreach (string key in dic_all.Keys)
            {
                file.WriteLine(key + "," + dic_all[key]);
            }
        }
    }
}
