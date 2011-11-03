using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace pepXMLReader
{
    class Program
    {
        static void Main(string[] args)
        {
            List <string> finalList = new List<string>();
            string currentDirectory = Directory.GetCurrentDirectory();
            string result = currentDirectory + "\\RESULT+3.csv";
            string[] filePaths = Directory.GetDirectories(currentDirectory);
            for (int i = 0; i < filePaths.Length; i++)
            {
                //Console.WriteLine("starting folder: " + filePaths[i]);
                string naiveFolder = filePaths[i] + "\\mvh-naive";
                string basoFolder = filePaths[i] + "\\mvh-baso";
                if (Directory.Exists(naiveFolder))
                {
                    string[] files = Directory.GetFiles(naiveFolder, "*.pepXML");
                    for (int k = 0; k < files.Length; k++)
                    {
                        string output = naiveFolder + "\\" + Path.GetFileNameWithoutExtension(files[k]) + "+3.csv";
                        string naivelist = files[k] + "," + "naive" + "," + Package.pepXMLReader(files[k], output);
                        finalList.Add(naivelist);
                    }
                }

                if (Directory.Exists(basoFolder))
                {
                    string[] files = Directory.GetFiles(basoFolder, "*.pepXML");
                    for (int k = 0; k < files.Length; k++)
                    {
                        string output = basoFolder + "\\" + Path.GetFileNameWithoutExtension(files[k]) + "+3.csv";
                        string basolist = files[k] + "," + "baso" + "," + Package.pepXMLReader(files[k], output);
                        finalList.Add(basolist);
                    }
                }
            }
            TextWriter tw = new StreamWriter(result);
            tw.WriteLine("file,model,num");
            foreach (string ss in finalList)
                tw.WriteLine(ss);
            tw.Flush();
            tw.Close();
        }
    }
}
