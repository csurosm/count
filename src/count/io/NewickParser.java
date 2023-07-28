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


import java.io.BufferedReader;
import java.io.IOException;
import java.io.PushbackReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import count.ds.IndexedTree;
import count.ds.Phylogeny;
import java.io.PrintStream;


/**
 * Newick-format phylogeny parsing.
 * 
 * Implementing classes fill exploit the hooks {@link #startChildren()}, 
 * {@link #nextChild() }, {@link #endChildren() } and {@link #endTree() },
 * {@link #setLength(double) } and {@link #setName(java.lang.String) }, 
 * which are called by {@link #readTree() }.
 *  The grammar used here is the following.
 * [Corresponds to the Newick format used by Phylip, DRAWTREE etc., with the addition of the '#' style comments,
 * see Joe Felsenstein;s <a href="http://evolution.genetics.washington.edu/phylip/newick_doc.html">specification</a>].
 *
 * The original specification allowed internal taxon names only after the subtree,
 * here we accept the name also before.
 *
 * Terminals:
 * <ul>
 * <li>SEMICOLON <code>;</code></li>
 * <li>LPAREN <code>(</code></li>
 * <li>RPAREN <code>)</code></li>
 * <li>COMMA <code>,</code></li>
 * <li>COLON <code>:</code></li>
 * <li>SEMICOLON <code>;</code></li>
 * <li>LBRACKET <code>[</code></li>
 * <li>RBRACKET <code>]</code></li>
 * <li>QUOTE <code>'</code></li>
 * <li>DBLQUOTE <code>"</code></li>
 * <li>NUMBER  IEEE floating point value  (Inf, -Inf, NaN are ok)</li>
 * <li>ALPHANUM</li>
 * </ul>
 * Grammar:
 * <pre>
 *	&lt;Tree&gt;           ::= &lt;Node&gt; SEMICOLON
 *	&lt;Node&gt;           ::= &lt;Leaf&gt;|&lt;Internal&gt;
 *	&lt;Internal&gt;       ::= &lt;Name&gt; LPAREN &lt;Nodelist&gt; RPAREN &lt;Edge length&gt;
 *						| LPAREN &lt;NodeList&gt; RPAREN &lt;Edge Length&gt;
 *						| LPAREN &lt;NodeList&gt; RPAREN &lt;Name&gt; &lt;Edge Length&gt;
 *	&lt;Leaf&gt;           ::= &lt;Name&gt; &lt;Edge length&gt;|&lt;Edge length&gt;
 *	&lt;Nodelist&gt;       ::= &lt;Node&gt;|(&lt;Node&gt; COMMA &lt;Nodelist&gt;)
 *	&lt;Edge length&gt;    ::= |&lt;Nonzero edge&gt;
 *	&lt;Nonzero edge&gt;   ::= COLON NUMBER
 *	&lt;Name&gt;           ::= &lt;quoted or unquoted name&gt;
 * </pre>
 *
 * Whitespaces (newline, tab, space) are allowed anywhere
 * between tokens. Remarks are allowed where whitespace is allowed,
 * they either start with '#' and continue until the end of line or are
 * enclosed in brackets.
 * 
 * 
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 * @since 2001
 */
public class NewickParser
{
    /**
     * Creates a new parser instance. 
     * 
     * @param input input source
     */
    protected NewickParser(Reader input)
    {
        this(new PushbackReader(input), 
             false, // no nesting of Newick comments
             true,  // hashmark comments are recognized
             true); // unescaped single quotes are recognized in taxon name
    }

    /**
     * Full instantiation of a parser with a given source.
     * 
     * @param input input source
     * @param nested_comments_allowed whether NEwick comments can be nested
     * @param hashmark_comments_allowed whether <code>#</code> can start a comment
     * @param relaxed_name_parsing whether unescaped single quotes are ok in the name
     */
    protected NewickParser(PushbackReader input, boolean nested_comments_allowed, boolean hashmark_comments_allowed, boolean relaxed_name_parsing)
    {    
        this.input = input;
        this.nested_comments_allowed = nested_comments_allowed;
        this.hashmark_comments_allowed = hashmark_comments_allowed;
        this.relaxed_name_parsing = relaxed_name_parsing;
    }
    
    private NewickParser(){this(null,false,true,true);}
    
    private final boolean nested_comments_allowed;
    private final boolean hashmark_comments_allowed;
    /**
     * Whether finicky parser should complain about improper quote usage (unescaped quotes, e.g., 'Baker's yeast') and 
     * messy node naming. 
     */
    private final boolean relaxed_name_parsing;
    private final PushbackReader input;
    
    
    // ---------------------------------------------------------
    // ---- Lexical analysis: terminal symbols.
    // ---------------------------------------------------------
    
    /** Newick format terminal: quote
     */
    public static final char QUOTE='\'';
    /** Newick format terminal: left parenthesis (starts a new set of descendants)
     */
    public static final char LPAREN='('; 
    /** Newick format terminal: right parenthesis (after last child)
     */
    public static final char RPAREN=')';
    /** Newick format terminal: comma (separates child subtrees)
     */
    public static final char COMMA=','; 
    /** Newick format terminal: semicolon (terminates the tree)
     */
    public static final char SEMICOLON=';';
    /** Newick format terminal: left bracket (starts a comment)
     */
    public static final char LBRACKET='[';
    /** Newick format terminal: left bracket (ends a comment)
     */
    public static final char RBRACKET=']';
    /** Newick format terminal: double quote 
     */
    public static final char DBLQUOTE='"';
    /** Newick format terminal: hashmark
     */
    public static final char HASHMARK='#';
    /** Newick format terminal: backslash
     */
    public static final char BACKSLASH='\\';
    /** Newick format terminal: colon (starts length)
     */
    public static final char COLON=':';
    /** Newick format terminal: underscore (translates into whitespace)
     */
    public static final char UNDERSCORE='_';
    /**
     * Characters that need to be enclosed in double quotes according to the specification
     */
    private static final String NEED_QUOTE_FOR= ""+QUOTE+LPAREN+RPAREN+COMMA+SEMICOLON+COLON+LBRACKET+RBRACKET;
    
    public static String printTree(Phylogeny tree)
    {
    	return printTree(tree
    			, false // quotes only if necessary
    			, true // show edge lengths
    			, false // no line breaks
    			, false // do not show identifiers
    			);
    }
    
    public static String printTree(IndexedTree tree)
    {
    	if (tree instanceof Phylogeny) return printTree((Phylogeny)tree);
    	else return printTree(new Phylogeny(tree));
    }
    
    /**
     * One-line tree description in Newick format with comments at some nodes.   
     * 
     * @param tree
     * @param node_comments
     * @return
     */
    public static String printTree(IndexedTree tree, String[] node_comments)
    {
    	Phylogeny phylo = 
    			tree instanceof Phylogeny
    			?(Phylogeny)tree
    			:new Phylogeny(tree);
    	boolean quote_name = false;
    	boolean show_length = true;
    	boolean line_breaks = false;
    	 String print_tree 
         = printTree(null, phylo.getRootNode(), quote_name, show_length, line_breaks, node_comments)
         	.append(';').toString();
    	 return print_tree;
    }
    
    /**
     * Newick-formatted String for a phylogeny.
     * 
     * @param tree (possibly empty, with null root)
     * @param always_quote_name whether single quotes are placed around node names alwats, or only when necessary
     * @param show_edge_lengths whether edge lengths should be printed
     * @param line_breaks_after_each_node whether each node should be in a new line (or, rather, return a one-line string)
     * @param show_node_identifiers whether the output should include then node indexes (in comments)
     * @return a properly formatted String by the specification
     */
    public static String printTree(
            Phylogeny tree,
            boolean always_quote_name,
            boolean show_edge_lengths,
            boolean line_breaks_after_each_node,
            boolean show_node_identifiers)
    {
        Phylogeny.Node root = tree.getRootNode();
        String[] node_comments = null;
        if (show_node_identifiers)
        {
        	node_comments = new String[tree.getNumNodes()];
        	for (int node=0; node<node_comments.length; node++)
        	{
        		node_comments[node] = tree.getNode(node).getNodeIdentifier();
        	}
        }
        String print_tree 
                = printTree(null, root, always_quote_name, show_edge_lengths, line_breaks_after_each_node, node_comments)
                .append(";").toString();
        return print_tree;
    }

    /** 
     * The recursive printing procedure for printing the subtree (without closing semicolon).
     * 
     * @param sb may be null
     * @param node root of subtree
     * @param always_quote_name
     * @param show_edge_lengths
     * @param line_breaks_after_each_node
     * @param node_comments
     * @return sb appended with printing of this subtree
     */
    private static StringBuilder printTree(
            StringBuilder sb,
            Phylogeny.Node node,
            boolean always_quote_name,
            boolean show_edge_lengths,
            boolean line_breaks_after_each_node,
            String[] node_comments)
    {
//        final String preNode = "< ";
//        final String inNode = "= ";
//        final String postNode = " > ";
        final String sepNode = line_breaks_after_each_node?"\n":" ";
        if (sb==null) sb = new StringBuilder();
        
        if (node == null) return sb;
        // previsit node info if requested
        int num_children = node.getNumChildren();
        // parenthesized list of children
        if (num_children>0)
        {
            sb.append(LPAREN);
            for (int ci=0; ci<num_children; ci++)
            {
                if (ci>0)
                {
                    sb.append(COMMA)
                    .append(sepNode);
                }
//                // infix order visit 
//                if (show_node_identifiers)
//                {
//                    sb.append(LBRACKET)
//                    .append(inNode)
//                    .append(node.getNodeIdentifier())
//                    .append(RBRACKET);
//                }
                sb = printTree(sb, node.getChild(ci), always_quote_name,show_edge_lengths,line_breaks_after_each_node,node_comments);
            }        
            sb.append(RPAREN);
        }
        // postfix visit, before name and edge length
        if (node_comments != null)
        {
        	String com = node_comments[node.getIndex()];
        	if (com != null)
        	{
	            sb.append(LBRACKET)
	            .append(com)
	            .append(RBRACKET);
        	}
        }
        
        if (node.isLeaf() || node.getName() != null)
        {
            if (!node.isLeaf()) sb.append(' ');
            sb.append(formatName(node.getName(),always_quote_name));
        }
        
        if (show_edge_lengths && !node.isRoot())
        {
            sb.append(COLON)
            .append(toIEEEFormat(node.getLength()));
        }

        // postvisit
//        if (show_node_identifiers)
//        {
//            sb.append(LBRACKET)
//            .append(postNode)
//            .append(node.getNodeIdentifier())
//            .append(RBRACKET);
//        }

        return sb;
    }
    

    /**
     * Inserts quotes into/around a string following Newick specification.
     * 
     * Quotes are added if a taxon name contains a lexical token (parentheses, brackets, quotes) or white space.
     * Single quote in the name ("baker's yeast") are encoded by doubling the single quote: <code>'baker''s yeast'</code>.
     * 
     * @param name taxon name (null is OK, empty String "" will be returned)
     * @param always_quote whether forcing quotes in the returned name (or only if necessary)
     * @return properly formatted string to write in Newick file as node name. 
     */
    public static String formatName(String name, boolean always_quote)
    {
        final String EMPTY_QUOTED_NAME = ""+QUOTE+QUOTE;
        if (name==null) return (always_quote?EMPTY_QUOTED_NAME:"");

        char[] cname=name.toCharArray();
        boolean need_quote=always_quote;
        for (int i=0; !need_quote && i<cname.length; i++)
        {
            char c = cname[i];
            need_quote = NEED_QUOTE_FOR.indexOf(c)!=-1 || Character.isWhitespace(c);
        }
        if (need_quote)
        {
            StringBuilder sb=new StringBuilder();
            sb.append(QUOTE);
            // replace single ' in the name with '' in the output
            int from_idx=0;
            int to_idx;
            do 
            {
                to_idx=name.indexOf(QUOTE,from_idx);
                if (to_idx==-1)
                    sb.append(cname,from_idx,cname.length-from_idx);
                else
                {
                    sb.append(cname,from_idx,to_idx-from_idx+1);
                    sb.append(QUOTE);
                    from_idx=to_idx+1;
                }
            } while(to_idx != -1 && from_idx<cname.length);
            sb.append(QUOTE);
            return sb.toString();
        } else
            return name;
    }
        
    /**
     * Default rounding precision in {@link #toIEEEFormat(double) }
     */
    private static final int DECIMALS=5;
    private static final double LOG10=Math.log(10.);
    /**
     * Small length under which a "negligible" value is displayed.
     * 
     * {@link #SHORT_BRANCH} should not be larger than 0.1^{@link #DECIMALS}.
     */
    private static final double SHORT_BRANCH = 1e-8;
    /**
     * 
     * Multiplier applied to {@link #SHORT_BRANCH} to display even shorter distances. 
     */
    private static final double PRACTICALLY_ZERO_BRANCH_MULTIPLIER = 0.099; 
  
    /**
     * Powers of ten up to at least 1/{@link #SHORT_BRANCH}.
     */
    private static final double[] POW10 = {1.,10.,100.,1000.,10000.,100000.,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22,1e23,1e24,1e25,1e26,1e27,1e28,1e29,1e30,1e31};

    /**
     * A compact display format including extended numerical values.
     * 
     * @param d the value that needs to be displayed
     * @return a formatted string with default precision 
     */
    public static String toIEEEFormat(double d)
    {
        return toIEEEFormat(d, DECIMALS, SHORT_BRANCH);
    }
    
    /**
     * A compact display format for extended numerical values.
     * 
     * NaN, positive and negative infinity are recognized. The value
     * is shown with a precision defined by the number of decimal digits 
     * wanted after the decimal point. Very short edge lengths are replaced with a 
     * fixed, tiny positive.
     * 
     * @param d value to be displayed
     * @param decimals preferred precision for large values
     * @param too_short length cutoff for very short edges
     * @return short character string for edge length display and I/O
     */
    public static String toIEEEFormat(double d, int decimals, double too_short)
    {
        assert (decimals>=0 && decimals<POW10.length);
        
        if (Double.isNaN(d)) return "NaN";
        else if (d==Double.POSITIVE_INFINITY) return "Inf";
        else if (d==Double.NEGATIVE_INFINITY) return "-Inf";
        else if (d==0.) return "0";
        else
        {
            int magnitude=(int) (Math.log(d)/LOG10+0.01);
            String retval;
            
            if (magnitude>=-(decimals+1))
            {
                // rounding by the precision of decimals
                double rounding_factor = POW10[decimals];
                double r=((int)(d*rounding_factor+0.5))/rounding_factor;
                retval= Double.toString(r);
            } else if (d<too_short) {
                retval= Double.toString(PRACTICALLY_ZERO_BRANCH_MULTIPLIER*too_short);//(short_branch*0.5)+"";
            } else
            {
                // keep only one digit after decimal point: "0.0_0x" or "xe-yy" returned 
                double m = POW10[-magnitude-1];
                double r=((int)(d*m)+0.5)/m;
          
                retval= Double.toString(r);
            }
            //System.out.println("#**TN.tIEEE "+d+"\tmag "+magnitude+"\tr "+retval);
            return retval;
            //return d+"";
        }
    }    
    
    
    // --------------------------------------------------------
    // ----- Parsing 
    // --------------------------------------------------------

    /**
     * Parsing state during reading.
     */
    private static enum ParsingState
    {
        BEFORE_NODE, 
        WITHIN_NODE,
        AFTER_NODE,
        PARSE_END;
    }
    
    /**
     * Parsing state: current node.
     */
    private Phylogeny.Node current_node;
    
    /**
     * Parsing state: current tree.
     */
    private Phylogeny current_tree;
    
    /**
     * Currently parsed node. 
     * 
     * @return 
     */
    protected Phylogeny.Node getCurrentNode(){ return current_node;}
    
    public static Phylogeny readTree(Reader R) throws ParseException, IOException
    {
        NewickParser parser = new NewickParser(R);
        Phylogeny readTree = parser.readTree();
//        System.out.println("#**NP.rT leaves "+Arrays.toString(readTree.getLeafNames()));
        return readTree;
    }
    
    public static Phylogeny[] readAllTrees(Reader R) throws ParseException, IOException
    {
        NewickParser parser = new NewickParser(R);
    	
    	List<Phylogeny> input_trees = new ArrayList<>();
    	int c;
    	while ((c=parser.skipBlanksAndComments())!=-1)
    	{
    		parser.input.unread(c);
    		Phylogeny P = parser.readTree();
    		input_trees.add(P);
    	}
    	return input_trees.toArray(new Phylogeny[0]);
    	
    }
    /**
     * Lexical parsing for Newick format. 
     * 
     * The codeuses the hooks {@link #startChildren()}, {@link #nextChild() }, {@link #endChildren() } and {@link #endTree() },
     * {@link #setLength(double) } and {@link #setName(java.lang.String) } for actions in different parsing states. 
     * @return the parsed phylogeny (until the next semicolon in the input)
     * @throws count.io.NewickParser.ParseException if input format does not conform
     * @throws IOException if I/O problem with file (reading or access) 
     */
    public Phylogeny readTree() throws ParseException, IOException
    {
        int current_level=0; // root level
        int c;
        ParsingState parsing_state=ParsingState.BEFORE_NODE;
        
        current_tree = new Phylogeny();
        current_node = current_tree.getRootNode();
        
        do
        {
            c=skipBlanksAndComments();
            
            //
            // --------------------------- LPAREN
            //
            if (c==LPAREN)
            {
                if (parsing_state == ParsingState.BEFORE_NODE)
                {
                    startChildren();
                    ++current_level;
                    // parsing_state=PARSE_BEFORE_NODE;
                } else
                    throw new ParseException(1, "Cannot have ``"+LPAREN+"'' here.");
            } else 
            //
            // --------------------------- COMMA
            //
            if (c==COMMA)
            {
                if (parsing_state == ParsingState.AFTER_NODE 
                        || parsing_state == ParsingState.WITHIN_NODE 
                        || parsing_state == ParsingState.BEFORE_NODE) 
                {
                    if (current_level==0)
                        throw new ParseException(2, "Cannot have ``"+COMMA+"'' at root level.");
                    nextChild();
                    parsing_state = ParsingState.BEFORE_NODE;
                } else
                    throw new ParseException(3, "Cannot have ``"+COMMA+"'' here.");
            } else
            //
            // --------------------------- RPAREN
            //
            if (c==RPAREN)
            {
                if (parsing_state == ParsingState.AFTER_NODE || parsing_state == ParsingState.WITHIN_NODE || parsing_state == ParsingState.BEFORE_NODE)
                {
                    if (current_level==0)
                        throw new ParseException(4, "Too many ``"+RPAREN+"''.");
                    --current_level;
                    endChildren();
                    parsing_state = ParsingState.WITHIN_NODE;
                } else
                    throw new ParseException(5, "Cannot have ``"+RPAREN+"'' here.");
            } else        
            //
            // --------------------------- COLON
            //
            if (c==COLON)
            {
                if (parsing_state == ParsingState.BEFORE_NODE || parsing_state == ParsingState.WITHIN_NODE)
                {
                    double d=parseEdgeLength();
                    setLength(d);
                    parsing_state=ParsingState.AFTER_NODE;
                } else
                    throw new ParseException(7,"Cannot have ``"+COLON+"'' here.");
            } else
            //
            // --------------------------- SEMICOLON
            //
            if (c==SEMICOLON)
            {
                if (parsing_state == ParsingState.AFTER_NODE || parsing_state == ParsingState.WITHIN_NODE || parsing_state == ParsingState.BEFORE_NODE)
                {
                    if (current_level != 0)
                        throw new ParseException(8,"Found ``"+SEMICOLON+"'' too early.");
                    endTree();
                    parsing_state=ParsingState.PARSE_END;
                }
            } else 
            //
            // --------------------------- taxon name
            //
            if (c!=-1)
            {
                if (parsing_state == ParsingState.WITHIN_NODE || parsing_state == ParsingState.BEFORE_NODE)
                {
//                    if (!relaxed_name_parsing && current_node.getName() != null)
//                        throw new ParseException(9, "Cannot name a node twice.");
                    input.unread(c);
                    String s=parseName();
                    setName(s);
                } else
                    throw new ParseException(10, "Cannot have node name here."); 
            }
            
        } while (c != -1 && parsing_state != ParsingState.PARSE_END);        
        
        
        if (parsing_state != ParsingState.PARSE_END)
            throw new ParseException(11, "Missing semicolon at the end");
        

        assert (c==SEMICOLON);
        
        return current_tree;
    }
    
    
    
    /**
     * Called by {@link #readTree() }  when starting to visit the current node's children (opening parenthesis).
     */
    protected void startChildren() 
    {
        current_node = current_node.newChild();
    }
    
    protected void nextChild() 
    {
        current_node = current_node.getParent().newChild();
    }

    protected void endChildren() 
    {
        current_node = current_node.getParent();
    }

    protected void setLength(double d) 
    {
        current_node.setLength(d);
        current_tree.hasLength(true);
    }
    
    protected void setName(String s) throws NewickParser.ParseException
    {
        if (current_node.getName()==null)
        {
            current_node.setName(s);
        } else if (relaxed_name_parsing)
        {
            current_node.setName(current_node.getName().concat(""+s));
        } else // 
        {
            assert (!relaxed_name_parsing && current_node.getName()!=null);
            throw new NewickParser.ParseException(9, "Cannot name a node twice.");
        } 
        
//        System.out.println("#**NP.setName "+s+"\t"+current_node);
    }

    protected void endTree()
    {
        // reset root to null if nothing there: a Newick file could describe an empty tree (";"). 
        assert (current_node == current_tree.getRootNode());
        Phylogeny.Node root = current_node;
        if (root.isLeaf() && root.getName()==null)
        {
            current_tree = new Phylogeny();
            current_node = null;
        }
    }
    
    
    /**
     * Set by {@link #skipBlanksAndComments()}
     */
    private StringBuilder bracketed_comment = null;
    
    /**
     * Skips white space and comments in the input. 
     * Reading position advances to the next informative character.
     * 
     * @return first non-comment, non-whitespace character
     * @throws IOException if reading fails
     */
    private int skipBlanksAndComments () throws IOException
    {
        int c;

        bracketed_comment = new StringBuilder();
        
        boolean parsed_a_comment;
        do
        {
        	parsed_a_comment=false; // reset at each iteration, so that more than 1 consecutive comments are skipped here 
            do{c=input.read();} while(c!=-1 && Character.isWhitespace((char)c)); // skip leading blanks
            
            if (c==LBRACKET)
            {
                parsed_a_comment=true; 
                int nesting_level=1;
                do
                {
                    c=input.read();
                    if (c==(int)RBRACKET) nesting_level--;
                    else if (nested_comments_allowed && c==(int)LBRACKET)
                        nesting_level++;
                    else if (c!=-1)// including, normally, left bracket that stays in the comment
                    {
                    	bracketed_comment.append((char)c);
                    }
                } while (nesting_level != 0 && c != -1);
            } else if (hashmark_comments_allowed && c==HASHMARK)
            	do {c=input.read();} while (c!=-1 && c!='\n' && c!='\r');
        } while(parsed_a_comment && c!=-1);
        return c;
    }    
    
    /**
     * Buffer for parsing node names and edge lengths. Initialized for 256 chars and extended when necessary.
     */
    private char[] buffer=new char[256];
    /**
     * The number of filled positions in the buffer. 
     */
    private int buffer_length=0;

    private void resetBuffer(){ buffer_length=0;}
    
    /**
     * Test if buffer array needs to be extended.
     */
    private void checkBuffer()
    {
        if (buffer_length == buffer.length)
        {
            int new_capacity=2*buffer.length;
            char[] new_buffer=new char[new_capacity];
            System.arraycopy(buffer,0,new_buffer,0,buffer.length);
            buffer=new_buffer;
        }
    }

    /**
     * Adds one character in the buffer; extended if necessary.
     * @param c character to be added at the end
     */
    private void addToBuffer(char c)
    {
        checkBuffer();
        buffer[buffer_length++]=c;
    }    
    
    /**
     * Parses edge length. 
     * In addition to usual numerical values, accept <code>NaN</code>, <code>+Inf</code>, <code>-Inf</code> and <code>Inf</code>.  
     * 
     * @return edge length; maybe NaN, positive or negative infinity
     * @throws IOException
     * @throws count.io.NewickParser.ParseException 
     */
    private double parseEdgeLength() throws IOException, ParseException
    {
        int c=skipBlanksAndComments();

        resetBuffer();

        while (c!=-1 && !Character.isWhitespace((char)c) && c!=COMMA && c!=RPAREN && c!=SEMICOLON)
        {
            addToBuffer((char)c);
            c=input.read();
        }

        double retval=1.0;

        if (buffer_length == 0) retval=0.;
        if (buffer_length >= 3)
        {
            if (buffer[0]=='N' && buffer[1]=='a' && buffer[2]=='N' && buffer_length==3)
                retval=Double.NaN;
            if (buffer_length==4 && buffer[1]=='I' && buffer[2]=='n' && buffer[3]=='f')
            { 
                if (buffer[0]=='-')
                    retval=Double.NEGATIVE_INFINITY;
                else if (buffer[0]=='+')
                    retval=Double.POSITIVE_INFINITY;
            }
            if (buffer[0]=='I' && buffer[1]=='n' && buffer[2]=='f' && buffer_length==3)
                retval=Double.POSITIVE_INFINITY;
        }

        if (retval==1.0) // stayed the same --- no special values were seen
        {
            try
            {
                retval=Double.parseDouble(new String(buffer,0,buffer_length));
            } catch (NumberFormatException e)
            {
                throw new ParseException(99,"Cannot parse edge length: "+e.toString());
            }
        }
        if (c!=-1) input.unread(c);

        return retval;
    }        

    /**
     * Reads in a taxon name according to Newick rules. 
     * If relaxed name parsing is allowed ({@link #relaxed_name_parsing}), then unescaped single quotes are 
     * accepted (e.g., <q>Brewer's yeast</q>). 
     * 
     * @return a String containing the taxon name 
     * @throws IOException if there is a problem while reading
     */
    private String parseName() throws IOException
    {
        resetBuffer();
        char quote=0;
        int c=input.read();
        if (c==QUOTE || c==DBLQUOTE)
        {
            quote=(char)c;
            c=input.read();
        }

        if (quote==0) while(c != -1)
        {
            // Unquoted labels may not contain blanks, parentheses, square brackets,
            // single_quotes, colons, semicolons, or commas.
            if (Character.isWhitespace((char)c) || NEED_QUOTE_FOR.indexOf(c)>-1)
                break;
            if (c==UNDERSCORE) // replaced with space according to specification
              c=' ';
            addToBuffer((char)c);
            c=input.read();
        }
        else while(c != -1) // quoted
        {
            if (c==quote)
            {
                // check whether next is also a quote
                c=input.read();
                if (c!=quote && !relaxed_name_parsing)
                { // we're done
                    break;
                }
                if (c==quote)
                {
                    addToBuffer((char)c);
                    input.read();
                } else
                {
                    // relaxed parsing: not an escaped quote but maybe just a mistake
                    if (c==-1
                        || Character.isWhitespace((char)c)
                        || c==LPAREN
                        || c==RPAREN
                        || c==COLON
                        || c==SEMICOLON
                        || c==COMMA)
                    { // definitely end of name
                        break;
                    }
                    // otherwise it was a mistake
                    addToBuffer(quote);
                    // no need to read(): it's already done
                }
            } else
            { // not a quote
                addToBuffer((char)c);
                c=input.read();
            }
        }
        if (c!=-1) input.unread(c);
        return new String(buffer,0,buffer_length);
    }        
    
    /**
     * Exception class for parsing syntax violations.
     */
    public static class ParseException extends IOException
    {
        protected ParseException(int error_id, String s)
        {
          super("Phylogeny parsing error (Newick format) "+error_id+":"+s);
        }
    }
    
    private void mainmain(String[] args) throws Exception
    {
    	
    	CommandLine cli = new CommandLine(args, getClass(), 0);
        Phylogeny tree = cli.getTree();
        PrintStream out = System.out;
        
        out.println(CommandLine.getStandardHeader(this.getClass()));
        out.println(CommandLine.getStandardRuntimeInfo(this.getClass(), args));
        
        boolean line_breaks = false;
        
		String filter_file = cli.getOptionValue(CommandLine.OPT_FILTER);
		if (filter_file != null)
    	{
    		List<String> leaves_kept = new ArrayList<>();
    		BufferedReader R = new BufferedReader(GeneralizedFileReader.guessReaderForInput(filter_file));
    		String line;
    		do
    		{
    			line = R.readLine();
    			if (line != null)
    			{
        			if (line.length()==0 || line.startsWith("#"))
        				continue;
        			String[] fields = line.split("\\t");
        			String name  = fields[0];
        			leaves_kept.add(name);
    			}
    		} while (line != null);
    		R.close();
    		tree.filterLeaves(leaves_kept.toArray(new String[0])); 
    	}
    	String relabel_file = cli.getOptionValue(CommandLine.OPT_RELABEL);
    	if (relabel_file != null)
    	{
    		Map<String,String> subst = new HashMap<>();
    		
    		BufferedReader R = new BufferedReader(GeneralizedFileReader.guessReaderForInput(relabel_file));
    		String line;
    		do
    		{
    			line = R.readLine();
    			if (line != null)
    			{
        			if (line.length()==0 || line.startsWith("#"))
        				continue;
        			String[] fields = line.split("\\t");
        			String old_name  = fields[0];
        			String new_name = fields[1];
        			subst.put(old_name, new_name);
    			}
    		} while (line != null);
    		R.close();

    		for (Phylogeny.Node L: tree.getLeaves())
        	{
        		String name = L.getName();
        		if (subst.containsKey(name))
        		{
        			L.setName(subst.get(name));
        			subst.remove(name);
        		} else
        		{
        			out.println("# no relabeling for "+L.getName());
        		}
        	}
//    		for (String name: subst.keySet())
//    		{
//    			out.println("# not in tree "+name);
//    		}
    	} // relabel
    	
    	if (filter_file==null && relabel_file==null)
        {
            for (Phylogeny.Node node: tree.getNodes())
            {
                out.println("#** node "+node);
            }
        	line_breaks = true;
        }
        out.println(printTree(tree, false, true, line_breaks, false));
    	
    }
    
    /**
     * Main entry for testing. 
     * 
     * @param args 1 or 2 arguments: file and comma-separated list of terminal taxa
     * @throws Exception whenever it feels like it
     */
    public static void main(String[] args) throws Exception
    {
    	NewickParser P = new NewickParser();
    	
//        if (args.length <1 || args.length> 2)
//        {
//            System.err.println("Call as $0 file [list]\n\twhere list is a comma-separated list of terminal taxa\n\t;Output is the tree spanned by those taxa.");
//            System.exit(0);
//        }
//        String tree_file = cli.getTree()
//        NewickParser P = new NewickParser(new java.io.FileReader(tree_file));
        P.mainmain(args);
    }

}
