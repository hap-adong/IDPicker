using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Data;
using pwiz.CLI.proteome;
using System.Xml;

//there is a bug reading the xml. in this file, I got 11549 decoy=0 sets, but actually, it has 10041 ids. 
namespace Pro2RNASeq
{
    class Program
    {
        static void Main(string[] args)
        {
            //get protein list from idpicker file .xml
            //recorded are the protein locus
            //remove reversed sequences (decoy = 1)
            XmlTextReader reader = new XmlTextReader("X:\\wangd5\\idpXML_FDR0.05\\myrimatch\\report_pep_2\\report_pep_2-assembled.idpXML");
            string locus = "";
            string decoy = "";
            List<string> locusList = new List<string>();
            while (reader.Read())
            {
                if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("protein"))
                {
                    locus = reader.GetAttribute("locus");
                    decoy = reader.GetAttribute("decoy");
                    if (decoy == "0") { locusList.Add(locus); Console.WriteLine(locus); }
                }
                if (reader.NodeType.Equals(XmlNodeType.Element) && reader.Name.Equals("peptideIndex")) break;
            }
            Console.WriteLine(locusList.Count);
            Console.WriteLine(Package.removeDuplicate(locusList).Count);

            /*
            string line;
            List<string> coloList = new List<string>();
            TextReader file2 = new StreamReader("X:\\wangd5\\colorectal.txt");
            while ((line = file2.ReadLine())!=null)
            {
                string[] str = line.Split(' ');
                coloList.Add(str[1]);
            }

            int mm = 0;
            foreach (string ss in coloList)
            { 
                if (locusList.Contains(ss)) mm++;
            }

            Console.WriteLine(mm);
            */
            
            
            
            
            
            
            
            
            
            //read rpkm file
            //put dictionary with:
            //key = protein locus; value = rpkm value. 
            Dictionary<string, string> rpkmDict = new Dictionary<string, string>();
            TextReader file = new StreamReader("X:\\wangd5\\s1_microarray95.rpkm");
            DataTable table = CSV.CsvParser.Parse(file, true);
            foreach (DataRow dr in table.Rows)
            {
                string longstr = dr[0].ToString();
                string[] str = longstr.Split(new char[] { ' ', '\t' });
                rpkmDict.Add(str[0], str[1]);
            }
            
            //map protein list to RNA_seq, get those mapped, output the rpkm values for R analysis
            //TextWriter mm_file = new StreamWriter("Z:\\wangd5\\temp\\mm_all.csv");
            int matched = 0;
            int unmatched = 0;
            foreach (var pair in rpkmDict)
            {
                if (locusList.Contains(pair.Key)) matched++;
                else unmatched++;
            }

            Console.WriteLine(matched + "||" + unmatched);
        }
    }
}
