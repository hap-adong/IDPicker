using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace NistLibraryReader
{
    class Program
    {
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

        public static string peptideBondParse (string peptide, int bond)
        {
            int expecedBond = peptide.Length - 1;
            string composition = "";
            if (expecedBond >= bond)
            {
                string NSeq = peptide.Substring(0, bond);
                string CSeq = peptide.Substring(bond, peptide.Length - bond);
                int totalNL = bond;
                int totalCL = peptide.Length - bond;

                //int NL = bond;
                //int CL = peptide.Length - bond;
                int nr = parseAAResidues(NSeq, 'R');
                int nk = parseAAResidues(NSeq, 'K');
                int nh = parseAAResidues(NSeq, 'H');
                int nd = parseAAResidues(NSeq, 'D');
                int ne = parseAAResidues(NSeq, 'E');
                int ns = parseAAResidues(NSeq, 'S');
                int nt = parseAAResidues(NSeq, 'T');
                int nq = parseAAResidues(NSeq, 'Q');
                int nn = parseAAResidues(NSeq, 'N');
                int nc = parseAAResidues(NSeq, 'C');
                int ng = parseAAResidues(NSeq, 'G');
                int np = parseAAResidues(NSeq, 'P');
                int na = parseAAResidues(NSeq, 'A');
                int nl = parseAAResidues(NSeq, 'L');
                int nm = parseAAResidues(NSeq, 'M');
                int nf = parseAAResidues(NSeq, 'F');
                int nw = parseAAResidues(NSeq, 'W');
                int ny = parseAAResidues(NSeq, 'Y');
                int nv = parseAAResidues(NSeq, 'V');

                int cr = parseAAResidues(CSeq, 'R');
                int ck = parseAAResidues(CSeq, 'K');
                int ch = parseAAResidues(CSeq, 'H');
                int cd = parseAAResidues(CSeq, 'D');
                int ce = parseAAResidues(CSeq, 'E');
                int cs = parseAAResidues(CSeq, 'S');
                int ct = parseAAResidues(CSeq, 'T');
                int cq = parseAAResidues(CSeq, 'Q');
                int cn = parseAAResidues(CSeq, 'N');
                int cc = parseAAResidues(CSeq, 'C');
                int cg = parseAAResidues(CSeq, 'G');
                int cp = parseAAResidues(CSeq, 'P');
                int ca = parseAAResidues(CSeq, 'A');
                int cl = parseAAResidues(CSeq, 'L');
                int cm = parseAAResidues(CSeq, 'M');
                int cf = parseAAResidues(CSeq, 'F');
                int cw = parseAAResidues(CSeq, 'W');
                int cy = parseAAResidues(CSeq, 'Y');
                int cv = parseAAResidues(CSeq, 'V');

                composition =  nr + "," + nk + "," + nh + "," + nd + "," + ne + "," + ns + "," +
                               nt + "," + nq + "," + nn + "," + nc + "," + ng + "," +
                               np + "," + na + "," + nl + "," + nm + "," + nf + "," +
                               nw + "," + ny + "," + nv + "," +
                               cr + "," + ck + "," + ch + "," + cd + "," + ce + "," + cs + "," +
                               ct + "," + cq + "," + cn + "," + cc + "," + cg + "," +
                               cp + "," + ca + "," + cl + "," + cm + "," + cf + "," +
                               cw + "," + cy + "," + cv;
                return composition;
            }
            return "";
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

        static void Main(string[] args)
        {
            string file = "Z:\\home\\dwang\\fragmentation\\NIST\\test\\charge3binary.csv";
            TextReader tr = new StreamReader(file);
            string line;
            int lineNumber = 0;
            List<string> list = new List<string>();
            List<string> evalList = new List<string>();
            int error = 0;
            //Dictionary<string, int> peptideDic = new Dictionary<string, int>();
            //string pre_peptide = "";
            //int pre_bond = 0;
            while ((line = tr.ReadLine()) != null)
            {
                lineNumber++;
                if (lineNumber == 1)
                    continue;
                string[] strs = line.Split(',');
                string peptide = strs[0];
                int bond = Convert.ToInt16(strs[1]);
                string NSeq = peptide.Substring(0, bond);
                string CSeq = peptide.Substring(bond, peptide.Length - bond);
                int nr = parseAAResidues(NSeq, 'R');
                int nk = parseAAResidues(NSeq, 'K');
                int nh = parseAAResidues(NSeq, 'H');
                int cr = parseAAResidues(CSeq, 'R');
                int ck = parseAAResidues(CSeq, 'K');
                int ch = parseAAResidues(CSeq, 'H');
                int totalNL = bond;
                int totalCL = peptide.Length - bond;
                string partialComposition = nr + "," + nk + "," + nh + "," + totalNL + "," + cr + "," + ck + "," + ch + "," + totalCL;
                //make label
                int y1intensity = Convert.ToInt16(strs[10]);
                int b1intensity = Convert.ToInt16(strs[11]);
                int y2intensity = Convert.ToInt16(strs[12]);
                int b2intensity = Convert.ToInt16(strs[13]);
                int y3intensity = Convert.ToInt16(strs[14]);
                int b3intensity = Convert.ToInt16(strs[15]);

                int label;
                bool b1y2 = false;
                bool b2y1 = false;
                if (b1intensity == 1 || y2intensity == 1) b1y2 = true;
                if (b2intensity == 1 || y1intensity == 1) b2y1 = true;
                
                if (b1y2 && !b2y1) label = 1;
                else if (b2y1 && !b1y2) label = 3;
                else if (b1y2 && b2y1) label = 2;
                else { label = 0; error++; }

                string composition = partialComposition + "," +  b1intensity + "," + b2intensity + "," + y1intensity + "," + y2intensity;

                if (label != 0)
                {
                    string final = composition + "," + label;
                    list.Add(final);
                }

                //for evaluation purposes:
                
                #region evaluation of logits
                
                //double logit_full = 1.3955 * nr + 1.1353 * nk + 1.3151 * nh + 0.3643 * nd + 0.4083 * ne + 0.3387 * ns + 0.4318 * nt + 0.5710 * nq + 0.5302 * nn + 1.1092 * nc + 0.2320 * ng + 0.5833 * np + 0.3368 * na + 0.4566 * nl + 0.6147 * nm + 0.5334 * nf + 0.6699 * nw + 0.5529 * ny + 0.4509 * nv -
                               //1.6576 * cr - 1.1609 * ck - 0.9583 * ch - 0.4717 * cd - 0.5540 * ce - 0.3934 * cs - 0.4631 * ct - 0.7360 * cq - 0.6134 * cn - 1.1476 * cc - 0.3483 * cg - 0.6572 * cp - 0.3537 * ca - 0.4927 * cl - 0.7008 * cm - 0.5854 * cf - 0.8251 * cw - 0.6350 * cy - 0.4399 * cv;
                double y1_logit = 0.1098112 * nr + 0.2085831 * nk + 0.1512109 * nh + 0.0460839 * totalNL - 0.3872417 * cr - 0.3684911 * ck - 0.1634741 * ch - 0.1693931 * totalCL + 1.2632997;
                double y2_logit = -0.6345364 * nr - 0.3365917 * nk - 0.4577882 * nh - 0.1492703 * totalNL + 0.7738133 * cr + 0.6036758 * ck + 0.5942542 * ch + 0.0701467 * totalCL + 0.0806280;
                double b1_logit = 0.0801432 * nr - 0.1088081 * nk - 0.1338220 * nh - 0.1413059 * totalNL - 0.3157957 * cr - 0.2708274 * ck - 0.3703136 * ch + 0.0157418 * totalCL + 1.2124699;
                double b2_logit = 0.8606449 * nr + 0.2763119 * nk + 0.4969152 * nh + 0.0685712 * totalNL - 1.3346995 * cr - 1.0977316 * ck - 1.0973677 * ch - 0.2028884 * totalCL + 1.9355980;
                double sd = nr * 0.9862 + nh * 0.8772 + nk * 0.7064 + totalNL * 0.4133 + cr * -1.1688 + ch * -0.3948 + ck * -0.6710 + totalCL * -0.4859;
                string evalFinal = b1intensity + "," + b2intensity + "," + y1intensity + "," + y2intensity + "," + b1_logit + "," + b2_logit + "," + y1_logit + "," + y2_logit + "," + sd;
                if (label != 0)
                    evalList.Add(evalFinal);
                 
                #endregion

            }
            tr.Close();

            Console.WriteLine(error);
            //string file_wr = "Z:\\home\\dwang\\fragmentation\\NIST\\test\\fullcomposition_binary.csv";
            //TextWriter tw = new StreamWriter(file_wr);
            //tw.WriteLine("nr,nk,nh,nd,ne,ns,nt,nq,nn,nc,ng,np,na,nl,nm,nf,nw,ny,nv,cr,ck,ch,cd,ce,cs,ct,cq,cn,cc,cg,cp,ca,cl,cm,cf,cw,cy,cv,b1,b2,y1,y2,label");
            //foreach (string str in list)
            //    tw.WriteLine(str);
            //tw.Flush();
            //tw.Close();

            string file_wr2 = "Z:\\home\\dwang\\fragmentation\\NIST\\test\\rocr.csv";
            TextWriter tw2 = new StreamWriter(file_wr2);
            tw2.WriteLine("b1,b2,y1,y2,lb1,lb2,ly1,ly2,sd");
            foreach (string str in evalList)
                tw2.WriteLine(str);
            tw2.Flush();
            tw2.Close();
        }
    }
}
