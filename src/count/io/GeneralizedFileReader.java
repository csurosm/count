/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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
package count.io;

import java.io.IOException;
import java.io.Reader;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

/**
 *
 * Common interface to reading files and URLs, possibly compressed.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public class GeneralizedFileReader extends java.io.BufferedReader
{
    
    public GeneralizedFileReader(String path_name) throws IOException
    {
        super(guessReaderForInput(path_name));
    }
    /**
     * Sets an input reader by parsing the file path: if it looks like an URL,
     * a URL connection is initiated; if it ends with <tt>gz</tt>, then
     * it is uncompressed on the fly.
     *
     * @param file_name URL or file path name
     * @return reader for the (possibly uncompressed file content)
     * @throws IOException if URL access fails, or the file cannot be opened for reading
     */
    public static Reader guessReaderForInput(String file_name) throws IOException
    {
        java.io.InputStream base ;
        if (file_name.startsWith("ftp:") || file_name.startsWith("http:"))
        {
            URL url = new URL(file_name);

            base = url.openStream();
        } else
        {
            base=new java.io.FileInputStream(file_name);
        }

        if (file_name.endsWith(".gz"))
            base = new GZIPInputStream(base);
        else if (file_name.endsWith("zip"))
            base = new ZipInputStream(base);
        return new java.io.InputStreamReader(base);
    }
}