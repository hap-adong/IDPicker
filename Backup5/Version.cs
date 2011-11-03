﻿//
// $Id: Version.cs 221 2010-11-04 21:22:23Z chambm $
//
// The contents of this file are subject to the Mozilla Public License
// Version 1.1 (the "License"); you may not use this file except in
// compliance with the License. You may obtain a copy of the License at
// http://www.mozilla.org/MPL/
//
// Software distributed under the License is distributed on an "AS IS"
// basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
// License for the specific language governing rights and limitations
// under the License.
//
// The Original Code is the IDPicker project.
//
// The Initial Developer of the Original Code is Matt Chambers.
//
// Copyright 2010 Vanderbilt University
//
// Contributor(s): Surendra Dasari
//

using System;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;

namespace IDPicker
{
    public static partial class Util
    {
        public static string Version { get { return GetAssemblyVersion(Assembly.GetExecutingAssembly().GetName()); } }
        public static DateTime LastModified { get { return GetAssemblyLastModified(Assembly.GetExecutingAssembly().GetName()); } }

        public static AssemblyName GetAssemblyByName (string assemblyName)
        {
            if (Assembly.GetCallingAssembly().GetName().FullName.Contains(assemblyName))
                return Assembly.GetCallingAssembly().GetName();

            foreach (AssemblyName a in Assembly.GetCallingAssembly().GetReferencedAssemblies())
            {
                if (a.FullName.Contains(assemblyName + ','))
                    return a;
            }
            return null;
        }

        public static string GetAssemblyVersion (AssemblyName assembly)
        {
            Match versionMatch = Regex.Match(assembly.ToString(), @"Version=([\d.]+)");
            return versionMatch.Groups[1].Success ? versionMatch.Groups[1].Value : "unknown";
        }

        public static DateTime GetAssemblyLastModified (AssemblyName assembly)
        {
            return File.GetLastWriteTime(Assembly.ReflectionOnlyLoad(assembly.FullName).Location);
        }
    }
}