package count.io;
/*
 * Copyright 2025 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Static methods for adding COG/arCOG/asCOG annotations.
 * 
 */
public class CogAnnotator {
	
	private static final String arcogs = "cogdata/arCOGs_names_Feb2022.txt";
	private static final String cogs = "cogdata/cog-24.def.tab";
	private static final String ascogs = "cogdata/asCOGs.2020-10.def.tab";
	
	/*
	 * Storing three Strings per family
	 */
	private static final int IDX_CATEGORY = 0;
	private static final int IDX_GENENAME = 1;
	private static final int IDX_ANNOTATION = 2;
	
	private static final String HDR_FAM = "OG name";
	private static final String HDR_CATEGORY = "Category";
	private static final String HDR_GENENAME = "Gene";
	private static final String HDR_ANNOTATION = "Description";
	
	
	
	private static final Map<String, String[]> og_annotations = new HashMap<>();
	private static boolean og_initialized = false;
	
	
	public static String[] getAvailableProperties() {
		if (!og_initialized) populateDataStructures() ;
		String[] getAvailableProperties = null;
		if (0<og_annotations.size()) {
			getAvailableProperties = new String[3];
			getAvailableProperties[IDX_CATEGORY] = HDR_CATEGORY;
			getAvailableProperties[IDX_GENENAME] = HDR_GENENAME;
			getAvailableProperties[IDX_ANNOTATION] = HDR_ANNOTATION;
		}
		return getAvailableProperties;
	}
	
	public static String[] getAnnotation(String og_name) {
		return og_annotations.get(og_name);
	}

	public static String[] getNullAnnotation() {
		String[] annot = new String[3];
		annot[IDX_CATEGORY] = "S";
		annot[IDX_GENENAME] = "";
		annot[IDX_ANNOTATION] = "Unknown";
		return annot;
	}
	
	private static boolean populateDataStructures() {
		IOException error1 = readAnnotations(arcogs,1,2,3);
		IOException error2 = readAnnotations(cogs,1,3,2);
		IOException error3 = readAnnotations(ascogs, 2,3,4);
		
		og_initialized = true;
		return error1==null && error2==null && error3==null;
	}

	private static IOException readAnnotations(String filename, int col_cat, int col_gene, int col_annot) {
		int col_fam = 0;
		
		try {
			URL url = ClassLoader.getSystemResource(filename);
			BufferedReader buf = new BufferedReader(new InputStreamReader(url.openStream()));
			String line = buf.readLine();
			while (line != null) {
				String[] fields = line.split("\t");
				if (col_annot<fields.length) {
					String fam = fields[col_fam];
					String[] annot = new String[3];
					annot[IDX_CATEGORY] = fields[col_cat];
					if (col_gene<fields.length)
						annot[IDX_GENENAME] = fields[col_gene];
					annot[IDX_ANNOTATION] = fields[col_annot];
					og_annotations.put(fam, annot);
				} else {
					//System.out.println("#**CogA.rA file "+filename+": skip line "+line);
				}
				line = buf.readLine();
			}
		} catch (IOException cannot_do_it) {
			return cannot_do_it;
		}
		return null;
	}
	

}
