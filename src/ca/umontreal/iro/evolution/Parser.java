package ca.umontreal.iro.evolution;

/**
 * Title:        Rooted tree reader class
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author       Miklos Csuros csuros AT iro.umontreal.ca
 * @version 1.0
 */

import java.io.IOException;
import java.io.Reader;
import java.io.PushbackReader;
import java.util.HashSet;
import java.util.Vector;

public class Parser {

  private Parser() {} // no instantiation

  // Newick format terminals
  public static final char QUOTE='\'';
  public static final char LPAREN='(';
  public static final char RPAREN=')';
  public static final char COMMA=',';
  public static final char SEMICOLON=';';
  public static final char LBRACKET='[';
  public static final char RBRACKET=']';
  public static final char DBLQUOTE='"';
  public static final char HASHMARK='#';
  public static final char BACKSLASH='\\';
  public static final char COLON=':';
  public static final char UNDERSCORE='_';

  public static final String NEED_QUOTE_FOR= ""+QUOTE+LPAREN+RPAREN+COMMA+SEMICOLON+COLON+LBRACKET+RBRACKET;

  /**
    *  The grammar used here is the following.
    * [Corresponds to the Newick format used by Phylip, DRAWTREE etc., with the addition of the '#' style comments,
    * see http://evolution.genetics.washington.edu/phylip/newick_doc.html].
    *
    * The original specification allowed internal taxon names only after the subtree,
    * here we accept the name before it also.
    *
    * Terminals:
    * SEMICOLON ;
    * LPAREN (
    * RPAREN )
    * COMMA ,
    * COLON :
    * SEMICOLON ;
    * LBRACKET [
    * RBRACKET ]
    * QUOTE '
    * DBLQUOTE "
    * NUMBER  IEEE floating point value  (Inf, -Inf, NaN are ok)
    * ALPHANUM
    *
    *	<Tree>           ::= <Node> SEMICOLON
    *	<Node>           ::= <Leaf>|<Internal>
    *	<Internal>       ::= <Name> LPAREN <Nodelist> RPAREN <Edge length>
    *						| LPAREN <NodeList> RPAREN <Edge Length>
    *						| LPAREN <NodeList> RPAREN <Name> <Edge Length>
    *	<Leaf>           ::= <Name> <Edge length>|<Edge length>
    *	<Nodelist>       ::= <Node>|(<Node> COMMA <Nodelist>)
    *	<Edge length>    ::= |<Nonzero edge>
    *	<Nonzero edge>   ::= COLON NUMBER
    *	<Name>           ::= <quoted or unquoted name>
    *
    * Whitespaces (newline, tab, space) are allowed anywhere
    * between tokens. Remarks are allowed where whitespace is allowed,
    * they either start with '#' and continue until the end of line or are
    * enclosed in brackets.

    */

  public static boolean NESTED_COMMENTS_ALLOWED=false;
  public static boolean HASHMARK_COMMENTS_ALLOWED=true;
  public static boolean RELAXED_NAME_PARSING=true;
  
  public static TreeNode readNewick(Reader input) throws IOException, ParseException 
  {
    return readNewick(input, NESTED_COMMENTS_ALLOWED, HASHMARK_COMMENTS_ALLOWED, RELAXED_NAME_PARSING);
  }
  
  public static TreeNode readNewick(
    Reader input,
    boolean nested_comments_allowed, // in PAUP, they are
    boolean hashmark_comments_allowed,
    boolean relaxed_name_parsing // tries to guess unquoted node names intelligently
    ) throws IOException, ParseException
  {
    TreeNode current_node=new TreeNode();
    TreeNode root=current_node;
    int c;
    int parsing_state=PARSE_BEFORE_NODE;

    PushbackReader R=new PushbackReader(input);
    
    do
    {
      c=skipBlanksAndComments(R,nested_comments_allowed,hashmark_comments_allowed);
      //
      // --------------------------- LPAREN
      //
      if (c==LPAREN)
      {
        if (parsing_state == PARSE_BEFORE_NODE)
        {
          current_node=current_node.newChild();
          // parsing_state=PARSE_BEFORE_NODE;
        } else
          throw new ParseException(1, "Cannot have ``"+LPAREN+"'' here.");
      } else
      //
      // --------------------------- COMMA
      //
      if (c==COMMA)
      {
        if (parsing_state == PARSE_AFTER_NODE || parsing_state == PARSE_WITHIN_NODE || parsing_state == PARSE_BEFORE_NODE)
        {
          if (current_node.isRoot())
            throw new ParseException(2, "Cannot have ``"+COMMA+"'' at root level.");
          current_node=current_node.mamma.newChild();
          parsing_state = PARSE_BEFORE_NODE;
        } else
          throw new ParseException(3, "Cannot have ``"+COMMA+"'' here.");
      } else
      //
      // --------------------------- RPAREN
      //
      if (c==RPAREN)
      {
        if (parsing_state == PARSE_AFTER_NODE || parsing_state == PARSE_WITHIN_NODE || parsing_state == PARSE_BEFORE_NODE)
        {
          if (current_node.isRoot())
            throw new ParseException(4, "Too many ``"+RPAREN+"''.");
          current_node=current_node.mamma;
          parsing_state = PARSE_WITHIN_NODE;
        } else
          throw new ParseException(5, "Cannot have ``"+RPAREN+"'' here.");
      } else
      //
      // --------------------------- COLON
      //
      if (c==COLON)
      {
        if (parsing_state == PARSE_BEFORE_NODE || parsing_state == PARSE_WITHIN_NODE)
        {
          double d=parseEdgeLength(R,nested_comments_allowed,hashmark_comments_allowed);
          current_node.setLength(d);
          parsing_state=PARSE_AFTER_NODE;
        } else
          throw new ParseException(7,"Cannot have ``"+COLON+"'' here.");
      } else
      //
      // --------------------------- SEMICOLON
      //
      if (c==SEMICOLON)
      {
        if (parsing_state == PARSE_AFTER_NODE || parsing_state == PARSE_WITHIN_NODE || parsing_state == PARSE_BEFORE_NODE)
        {
          if (!current_node.isRoot())
            throw new ParseException(8,"Found ``"+SEMICOLON+"'' too early.");
          parsing_state=PARSE_END;
        }
      } else if (c!=-1)
      {
      //
      // --------------------------- taxon name
      //
        if (parsing_state == PARSE_WITHIN_NODE || parsing_state == PARSE_BEFORE_NODE)
        {
          if (!relaxed_name_parsing && current_node.name != null)
            throw new ParseException(9, "Cannot name a node twice.");
          R.unread(c);
          String s=parseName(R,relaxed_name_parsing);
          if (current_node.name != null)
            current_node.name=current_node.name.concat(" && "+s);
          else
            current_node.name=s;
        } else
          throw new ParseException(10, "Cannot have node name here.");
      }
    } while (c != -1 && parsing_state != PARSE_END);
    
    if (parsing_state == PARSE_BEFORE_NODE && root.isLeaf())
        return null;
    
    // check if all leaves have a name
    TreeNode[] leaves = root.getTraversal().getLeaves();
    for (int leaf_idx=0; leaf_idx<leaves.length; leaf_idx++)
        if (leaves[leaf_idx].getName()==null)
            throw new ParseException(11, "Tree has unnamed terminal nodes (too many commas somewhere maybe?)");
    
    return root;
  }

  // parsing states
  private static final int PARSE_BEFORE_NODE=1;
  private static final int PARSE_WITHIN_NODE=2;
  private static final int PARSE_AFTER_NODE=3;
  private static final int PARSE_END=4;
  
  public static TreeNode[] readTrees
          (
            Reader input,
            boolean nested_comments_allowed, // in PAUP, they are
            boolean hashmark_comments_allowed,
            boolean relaxed_name_parsing // tries to guess unquoted node names intelligently
           ) throws IOException, ParseException
  {
    Vector<TreeNode> trees = new Vector<TreeNode>();
    while (true)
    {
        TreeNode R = readNewick(input, nested_comments_allowed, hashmark_comments_allowed, relaxed_name_parsing);
        if (R==null)
            break;
        trees.add(R);
    }
    return trees.toArray(new TreeNode[0]);  
  }

  public static TreeNode[] readTrees(Reader input) throws IOException, ParseException 
  {
    return readTrees(input, NESTED_COMMENTS_ALLOWED, HASHMARK_COMMENTS_ALLOWED, RELAXED_NAME_PARSING);
  }

  public static int nextNonWhitespace(PushbackReader input) throws IOException
  {
    return skipBlanksAndComments(input, false, false);
  }

  private static int skipBlanksAndComments (
    PushbackReader input,
    boolean nested_comments_allowed,
    boolean hashmark_comments_allowed) throws IOException
  {
  return skipBlanksAndComments(input, true, nested_comments_allowed, hashmark_comments_allowed);
  }
  private static int skipBlanksAndComments (
    PushbackReader input,
    boolean comments_allowed,
    boolean nested_comments_allowed,
    boolean hashmark_comments_allowed) throws IOException
  {
    int c;
    double d=1.0;
    boolean parsed_a_comment=false;
    do
    {
      do{c=input.read();} while(c!=-1 && Character.isWhitespace((char)c)); // skip leading blanks
      if (comments_allowed && c==LBRACKET)
      {
        parsed_a_comment=true;
        int nesting_level=1;
        do
        {
          c=input.read();
          if (c==RBRACKET) nesting_level--;
          if (nested_comments_allowed && c==LBRACKET)
            nesting_level++;
        } while (nesting_level != 0 && c != -1);
      } else if (hashmark_comments_allowed && c==HASHMARK)
        do {c=input.read();} while (c!=-1 && c!='\n' && c!='\r');
    } while(parsed_a_comment && c!=-1);

    return c;
  }

  private static int buffer_capacity=256;
  private static char[] buffer=new char[buffer_capacity];
  private static int buffer_length=0;

  private static void checkBuffer()
  {
    if (buffer_length == buffer_capacity)
    {
      int new_capacity=2*buffer_capacity;
      char[] new_buffer=new char[new_capacity];
      System.arraycopy(buffer,0,new_buffer,0,buffer_capacity);
      buffer=new_buffer;
      buffer_capacity=new_capacity;
    }
  }
  private static void addToBuffer(char c)
  {
    checkBuffer();
    buffer[buffer_length++]=c;
  }

  public static double parseDouble(PushbackReader input) throws IOException, ParseException
  {
    return parseEdgeLength(input,false,false,false,false);
  }

  private static double parseEdgeLength(
    PushbackReader input,
    boolean nested_comments_allowed,
    boolean hashmark_comments_allowed) throws IOException, ParseException
  {
    return parseEdgeLength(input,true, nested_comments_allowed, hashmark_comments_allowed, true);
  }

  private static double parseEdgeLength(
    PushbackReader input,
    boolean comments_allowed,
    boolean nested_comments_allowed,
    boolean hashmark_comments_allowed,
    boolean whitespace_allowed) throws IOException, ParseException
  {
    int c;
    if (whitespace_allowed)
      c=skipBlanksAndComments(input,comments_allowed,nested_comments_allowed,hashmark_comments_allowed);
    else
      c=input.read();

    buffer_length=0;
    while (c!=-1 && !Character.isWhitespace((char)c) && c!=COMMA && c!=RPAREN && c!=SEMICOLON)
    {
      addToBuffer((char)c);
      c=input.read();
    }

    double retval=1.0;

    if (buffer_length == 0)
      retval=0.;

    if (buffer_length >= 3)
    {
      if (buffer[0]=='N' && buffer[1]=='a' && buffer[2]=='N' && buffer_length==3)
        retval=Double.NaN;
      if (buffer_length==4 && buffer[1]=='I' && buffer[2]=='n' && buffer[3]=='f')
      { if (buffer[0]=='-')
          retval=Double.NEGATIVE_INFINITY;
        else if (buffer[0]=='+')
          retval=Double.POSITIVE_INFINITY;
      }
      if (buffer[0]=='I' && buffer[1]=='n' && buffer[2]=='f' && buffer_length==3)
        retval=Double.POSITIVE_INFINITY;
    }

    if (retval==1.0)
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

  public static String parseString(
    PushbackReader input, boolean relaxed_name_parsing) throws IOException
  {
    return parseName(input,relaxed_name_parsing);
  }

  private static String parseName(
    PushbackReader input,
    boolean relaxed_name_parsing) throws IOException
  {
    buffer_length=0;
    char quote=0;
    int c=input.read();
    if (c==QUOTE || c==DBLQUOTE)
    {
      quote=(char)c;
      c=input.read();
    }

    if (quote==0)
      while(c != -1)
      {
        // Unquoted labels may not contain blanks, parentheses, square brackets,
        // single_quotes, colons, semicolons, or commas.
        if (Character.isWhitespace((char)c)
            || NEED_QUOTE_FOR.indexOf(c)>-1)
            break;
        if (c==UNDERSCORE) // replaced with space
          c=' ';
        addToBuffer((char)c);
        c=input.read();
      }
    else // quoted
      while(c != -1)
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


  public static class ParseException extends Exception
  {
    public ParseException(int error_id, String s)
    {
      super("Parsing error "+error_id+":"+s);
    }
  }
  
  /**
   * Computes the subtree spanned by a set of terminal taxa
   */
  public static TreeNode spanTree(TreeNode root, String[] selected_taxon_name)
  {
      HashSet<String> selected_leaves = new HashSet<String>();
      for (int i=0; i<selected_taxon_name.length; i++)
          selected_leaves.add(selected_taxon_name[i]);
      TreeNode[] nodes = root.subtreeNodes();
      for (int i=0; i<nodes.length; i++)
          nodes[i].setId(i);
      TreeNode[] new_nodes = new TreeNode[nodes.length];
      for (int node_idx=0; node_idx<nodes.length; node_idx++)
      {
          TreeNode old_node = nodes[node_idx];
          if (old_node.isLeaf())
          {
              String name = old_node.getName();
              if (selected_leaves.contains(name))
              {
                  // okay, we need this guy
                  new_nodes[node_idx]=old_node;
              } else
              {
                  new_nodes[node_idx]=null;
              }
          } else
          {
              int num_children_spanned = 0;
              int num_children = old_node.getNumChildren();
              for (int cidx=0; cidx<num_children; cidx++)
              {
                  TreeNode Child = old_node.getChild(cidx);
                  if (new_nodes[Child.getId()]!=null)
                  {
                      num_children_spanned++;
                  }
              }
              if (num_children_spanned == 0)
              {
                  new_nodes[node_idx]=null;
              } else if (num_children_spanned == 1)
              {
                  // find that one child and move it up into this slot
                  for (int cidx=0; cidx<num_children; cidx++)
                  {
                      TreeNode Child = old_node.getChild(cidx);
                      if (new_nodes[Child.getId()]!=null)
                      {
                          TreeNode newChild = new_nodes[Child.getId()];
                          new_nodes[node_idx]=newChild;
                          newChild.setLength(newChild.getLength()+old_node.getLength());
                          break;
                      }
                  }
              } else 
              {
                  TreeNode N = new TreeNode();
                  for (int cidx=0; cidx<num_children; cidx++)
                  {
                      TreeNode Child = old_node.getChild(cidx);
                      if (new_nodes[Child.getId()]!=null)
                      {
                          TreeNode newChild = new_nodes[Child.getId()];
                          N.addChild(newChild);
                      }
                  }
                  N.setLength(old_node.getLength());
                  new_nodes[node_idx]=N;
                  N.setId(node_idx);
              }
          }
      } // for node_idx
      return new_nodes[root.getId()];
  }
  
  /** Attempts to select the subtere spanned by a list of organisms
   */
  public static void main(String[] args) throws java.io.IOException, ParseException
  {
      if (args.length != 2)
      {
          System.err.println("Call as $0 file list\n\twhere list is a comma-separated list of terminal taxa\n\t;Output is the tree spanned by those taxa.");
          System.exit(0);
      }
      String tree_file = new String(args[0]);
      TreeNode root = readNewick(new java.io.FileReader(tree_file));

      String taxon_list = new String(args[1]);
      String[] selected_taxon = taxon_list.split(",");
      
      TreeNode new_root = spanTree(root,selected_taxon);
      System.out.println(new_root.newickTree(false, false, false, false));
  }

  static String RCS_ID="$Id: Parser.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: Parser.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //
}