
package ca.umontreal.iro.evolution.malin.ui.count;

import java.util.Enumeration;
import java.util.Hashtable;

import java.awt.Font;
import java.awt.Frame;

import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;

import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;

import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;


/**
 *
 * @author csuros
 */
public class UsersGuide extends JDialog
{

    /** Creates a new instance of UsersGuide */
    public UsersGuide(Frame owner)
    {
        super(owner, "User's Guide");
        initDataStructures();
        initComponents();
        pack();

    }


    private DefaultMutableTreeNode root;
    private DefaultTreeModel tree_model;

    private static final Object ROOT_USER_OBJECT = new Object();

    private Hashtable<String, DefaultMutableTreeNode> guide_toc;

    private void initDataStructures()
    {
        root = new DefaultMutableTreeNode(ROOT_USER_OBJECT,true);

        tree_model = new DefaultTreeModel(root);

        guide_toc = new Hashtable<String, DefaultMutableTreeNode>();

        createNodes();
    }

    private JTree navigation;
    private JSplitPane content;

    private static final EntryPoint EMPTY_ENTRY_POINT = new EntryPoint("[dummy]",null);

    private void initComponents()
    {
        navigation = new JTree(tree_model);

        navigation.setRootVisible(false);
        navigation.setShowsRootHandles(true);
        navigation.setEditable(false);
        navigation.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        navigation.addTreeSelectionListener(new TreeSelectionListener()
        {
            public void valueChanged(TreeSelectionEvent E)
            {
                setEntryPoint();
            }
        });

        navigation.setFont(new Font("Serif", Font.BOLD,14));
        Enumeration all_nodes = root.preorderEnumeration();
        while(all_nodes.hasMoreElements())
        {
            DefaultMutableTreeNode N = (DefaultMutableTreeNode)all_nodes.nextElement();
            navigation.expandPath(new TreePath(N.getPath()));
        }

        content = new JSplitPane();
        content.setDividerLocation(300+content.getInsets().top);
        content.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        content.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
//        content.setBorder(null);
        content.setResizeWeight(0.0);
        content.setOneTouchExpandable(true);

        content.setLeftComponent(new JScrollPane(navigation));
        content.setRightComponent(EMPTY_ENTRY_POINT.displayed_pane);

        add(content);
    }


    private void createNodes()
    {
        DefaultMutableTreeNode current_chapter = null;
        DefaultMutableTreeNode current_section = null;
        current_chapter =   newNode(root, CHAPTER_OVERVIEW_T, null);
        current_section =   newNode(current_chapter, SECTION_INTRODUCTION_T, SECTION_INTRODUCTION );
        current_section =   newNode(current_chapter, SECTION_PARSIMONY_T, SECTION_PARSIMONY );
        current_section =   newNode(current_chapter, SECTION_PHYLOBD_T, null );
                            newNode(current_section, SUBSECTION_PHYLOBD_RATES_T, SUBSECTION_PHYLOBD_RATES);
                            newNode(current_section, SUBSECTION_XENOLOGS_INPARALOGS_T, SUBSECTION_XENOLOGS_INPARALOGS);
                            newNode(current_section, SUBSECTION_POSTERIORS_T, SUBSECTION_POSTERIORS);
                            newNode(current_section, SUBSECTION_ABSENT_FAMILIES_T, SUBSECTION_ABSENT_FAMILIES);
        current_chapter =   newNode(root, CHAPTER_USING_COUNT_T, null);
        current_section =   newNode(current_chapter, SECTION_DESIGN_CONCEPTS_T, null);
                            newNode(current_section, SUBSECTION_ERRORS_T, SUBSECTION_ERRORS);
                            newNode(current_section, SUBSECTION_SESSIONS_WORKAREA_T, SUBSECTION_SESSIONS_WORKAREA);
                            newNode(current_section, SUBSECTION_BROWSERS_T, SUBSECTION_BROWSERS);
                            newNode(current_section, SUBSECTION_TREE_DISPLAYS_T, SUBSECTION_TREE_DISPLAYS);
                            newNode(current_section, SUBSECTION_TABLE_COLUMNS_T, SUBSECTION_TABLE_COLUMNS);
                            newNode(current_section, SUBSECTION_COMMENTS_T, SUBSECTION_COMMENTS);
        current_section =   newNode(current_chapter, SECTION_SESSIONS_T, null);
                            newNode(current_section, SUBSECTION_PHYLOGENY_T, SUBSECTION_PHYLOGENY);
                            newNode(current_section, SUBSECTION_SAVING_T, SUBSECTION_SAVING);
                            newNode(current_section, SUBSECTION_TABLE_T, SUBSECTION_TABLE);
                            newNode(current_section, SUBSECTION_ANNOTATIONS_T, SUBSECTION_ANNOTATIONS);
                            newNode(current_section, SUBSECTION_FILTERING_T, SUBSECTION_FILTERING);
        current_section =   newNode(current_chapter, SECTION_RATES_T, SECTION_RATES);
                            newNode(current_section, SUBSECTION_RATE_FILE_T, SUBSECTION_RATE_FILE);
                            newNode(current_section, SUBSECTION_RATES_PANEL_T, SUBSECTION_RATES_PANEL);
                            newNode(current_section, SUBSECTION_RATE_OPTIMIZATION_T, SUBSECTION_RATE_OPTIMIZATION);
        current_section =   newNode(current_chapter, SECTION_ANALYSIS_T, SECTION_ANALYSIS);
                            newNode(current_section, SUBSECTION_ANALYSIS_PANEL_T, SUBSECTION_ANALYSIS_PANEL);
                            newNode(current_section, SUBSECTION_ANALYSIS_DOLLO_T, SUBSECTION_ANALYSIS_DOLLO);
                            newNode(current_section, SUBSECTION_ANALYSIS_WAGNER_T, SUBSECTION_ANALYSIS_WAGNER);
                            newNode(current_section, SUBSECTION_ANALYSIS_POSTERIORS_T, SUBSECTION_ANALYSIS_POSTERIORS);
                            newNode(current_section, SUBSECTION_ANALYSIS_PGL_T, SUBSECTION_ANALYSIS_PGL);
        current_chapter =   newNode(root, CHAPTER_REFERENCES_T, CHAPTER_REFERENCES);

    }

    private DefaultMutableTreeNode newNode(DefaultMutableTreeNode parent, String title, String body)
    {
        EntryPoint E = new EntryPoint(title, body);
        DefaultMutableTreeNode N = new DefaultMutableTreeNode(E);
        parent.add(N);
        guide_toc.put(title, N);

        return N;
    }

    /**
     * Selects the node with this title
     */
    public void openPage(String title)
    {
        DefaultMutableTreeNode N = guide_toc.get(title);
        if (N != null)
            navigation.setSelectionPath(new TreePath(N.getPath()));
    }

    /**
     * Sets the displayed entry point based on the tree node selecion
     */
    private void setEntryPoint()
    {
        TreePath selected = navigation.getSelectionPath();
        int div = content.getDividerLocation(); // so that the divider doesn't change

        if (selected == null)
        {
            content.setRightComponent(EMPTY_ENTRY_POINT.displayed_pane);
        } else
        {
            DefaultMutableTreeNode selected_node = (DefaultMutableTreeNode)selected.getLastPathComponent();
            EntryPoint EP = (EntryPoint) selected_node.getUserObject();
            content.setRightComponent(EP.displayed_pane);
        }

        content.setDividerLocation(div);
    }








   /* ... */

   /**
    * Class for storing title-body parts of
    * sections, subsections, etc
    */
   private static class EntryPoint
   {
       private EntryPoint(String title, String body)
       {
           this.title = title;
           JEditorPane text_pane = new JEditorPane();
           text_pane.setContentType("text/html");
           text_pane.setEditable(false);
           if (body != null)
               text_pane.setText("<h2>"+title+"</h2>"+body);
           else
               text_pane.setText("");
           displayed_pane = new JScrollPane(text_pane);

       }


       private String title;

       private JScrollPane displayed_pane;

       public String toString()
       {
           return title;
       }
   }


    public static final String CHAPTER_OVERVIEW_T = "Overview";

    public static final String SECTION_INTRODUCTION_T = "Introduction";
    public static final String SECTION_INTRODUCTION
            = "<p>Count is a " +
            "software package for the evolutionary analysis of homolog " +
            "family sizes, or other numerical census-type characters " +
            "along a phylogeny.</p>" +
            "<p>The " +
            "methods implemented in Count aim to facilitate the " +
            "evolutionary analysis of gene content evolution. The " +
            "principal data consist of the distribution of homolog family " +
            "sizes across multiple genomes. In particular, the data " +
            "correspond to a table &Phi;<sub><var>f</var>,<var>j</var></sub>: " +
            "<var>f</var>=1,..., <var>n</var>; " +
            "<var>j</var>=1,..., <var>m</var>, where &Phi;<sub><var>f</var>,<var>j</var></sub> " +
            "is the number of " +
            "homologous genes to family <var>f</var> that are found in genome <var>j</var>. " +
            "The vector &Phi;<sub><var>f</var></sub>=(&Phi;<sub><var>f</var>,<var>j</var></sub>: <var>j</var>=1,..., <var>m</var>) " +
            "is the so-called <em>phylogenetic profile</em> of family <var>f</var>. " +
            "<em>Ancestral reconstruction</em> is " +
            "the problem of inferring family sizes at inner nodes of a " +
            "given evolutionary tree over a subset of the genomes " +
            "<var>j</var>=1,... , <var>m</var>. " +
            "In a <em>parsimony</em> " +
            "approach, the phylogenetic profile is extended to inner " +
            "nodes by minimizing a penalty function over the implied size" +
            "changes on the tree edges. " +
            "A <em>likelihood</em> approach assumes an explicit " +
            "probabilistic model for phylogenetic profiles. Count " +
            "implements so-called phylogenetic birth-and-death models, in " +
            "which family size evolution on an edge is governed by a " +
            "linear birth-and-death model traditionally employed in the " +
            "contexts of queuing systems and population growth. After " +
            "optimizing the parameters of such a model on a data set, " +
            "ancestral reconstruction can be carried out by computing " +
            "posterior probabilities for the family sizes at inner nodes.</p>" +
            "<p>For further details about the algorithmic procedures and the " +
            "mathematical framework, please consult the " +
            "<a href=\"http://www.iro.umontreal.ca/~csuros/gene_content/Count/count-usage.pdf\">PDF version of " +
            "the User's Guide</a>.</p>";

    public static final String SECTION_PARSIMONY_T = "Parsimony";
    public static final String SECTION_PARSIMONY
            = "<p>Count implements a parsimony method known as asymmetric Wagner " +
            "parsimony. The method is " +
            "described in details elsewhere [Cs&#369;r&ouml;s 2008]. " +
            "The key idea is to " +
            "penalize the changes of gene family size differently in " +
            "cases of losses and gains. Specifically, a change from <var>x</var> " +
            "to <var>y</var> on an edge is penalized either by (<var>x</var>-<var>y</var>) " +
            "when <var>x</var>&ge;<var>y</var>, " +
            "or by <var>g</var>(<var>y</var>-<var>x</var>) " +
            "when <var>y</var>&gt;<var>x</var>. " +
            "The <var>g</var> parameter sets the " +
            "relative penalty of a gain vs. a loss. " +
            "In classic Wagner parsimony [Farris, 1970], " +
            "<var>g</var>=1, but if losses happen more often than " +
            "gains, then some <var>g</var>&gt;1 may be a more adequate choice.</p>";
    public static final String SECTION_PHYLOBD_T = "Phylogenetic birth-and-death model";
    public static final String SUBSECTION_PHYLOBD_RATES_T = "Rates";
    public static final String SUBSECTION_PHYLOBD_RATES
            = "<p>A <em>phylogenetic birth-and-death model</em> " +
            "assumes that a stochastic process " +
            "acts on each edge, determining the evolution of homolog " +
            "family size. The process on an edge is characterized by " +
            "three parameters, denoted by &kappa;, &mu;,and &lambda;. " +
            "A family of size <var>n</var> decreases by a rate " +
            "of <var>n</var>&mu; " +
            "and increases by a rate of (&kappa;+<var>n</var>&lambda;). " +
            "In the context of a homolog gene family, &mu; is the " +
            "individual gene loss rate " +
            "(uniform across members of the family), " +
            "&lambda; is the individual gene " +
            "duplication rate (uniform across " +
            "members of the family), " +
            "and &kappa; is the rate of gene " +
            "gain by any mechanism, including innovation " +
            "and lateral gene transfer. " +
            "The model and the associated " +
            "computational techniques are described in details elsewhere " +
            "[Cs&#369;r&ouml;z &amp; Mikl&oacute;s 2006, 2009a, 2009b].</p>" +
            "<p>In the most general model, the process " +
            "parameters (&kappa;, &mu;, &lambda;) differ across edges, and " +
            "depend on the gene family. Specifically, the linear " +
            "birth-and-death process on edge <var>e</var> for family <var>f</var> " +
            "has rate parameters " +
            "&kappa;=&kappa;<sub><var>e</var></sub>&kappa;'<sub><var>f</var></sub>, " +
            "&mu;=&mu;<sub><var>e</var></sub>&mu;'<sub><var>f</var></sub>, " +
            "&lambda;=&lambda;<sub><var>e</var></sub>&lambda;'<sub><var>f</var></sub>, " +
            "and runs for a duration of <var>t</var><sub><var>e</var></sub> <var>t</var>'<sub><var>f</var></sub>. " +
            "Edge length is <var>t</var><sub><var>e</var></sub>; " +
            "&kappa<sub><var>e</var></sub>, &mu;<sub><var>e</var></sub>, &lambda;<sub><var>e</var></sub>  " +
            "are lineage-specific average rate parameters; " +
            "<var>t</var>'<sub><var>f</var></sub>, " +
            "&kappa;'<sub><var>f</var></sub>, " +
            "&mu;'<sub><var>f</var></sub>, " +
            "&lambda;'<sub><var>f</var></sub>, " +
            "are family-specific rate factors. " +
            "These latter are assumed to be either constant, or" +
            "have a discretized Gamma distribution [Yang 1994]. " +
            "For gain and duplication rates, there is a possibility for " +
            "mixing in 0-rate categories: &kappa;'<sub><var>f</var></sub>=0 or " +
            "&lambda;'<sub><var>f</var></sub>=0 " +
            "with some prior probabilities " +
            "set during optimization. Then the " +
            "family-specific duplication and gain rate factors have the " +
            "discretized Gamma distribution in the non-zero-rate " +
            "categories. Model parameters are set by likelihood " +
            "optimization in Count. Given model parameters yield " +
            "exactly computable likelihoods and posterior probabilities " +
            "for ancestral gene content: see Cs&#369;r&ouml;s &amp; Mikl&oacute;s " +
            "[2009a, 2009b].</p>";
    public static final String SUBSECTION_XENOLOGS_INPARALOGS_T = "Xenologs and inparalogs";
    public static final String SUBSECTION_XENOLOGS_INPARALOGS
            = "<p>A key notion in the " +
            "likelihood computation is that of partitioning the family " +
            "members at a child node into <em>xenolog</em> " +
            "and <em>inparalog</em>. " +
            "The xenolog group consists of the members that have no ancestor at the " +
            "parent node, i.e., their ancestor appeared in a gain event " +
            "within the lineage leading to the child node. Each family " +
            "member at the ancestor node has a corresponding inparalog " +
            "group formed by its descendants at the child node.<p>";

    public static final String SUBSECTION_POSTERIORS_T = "Posteriors";
    public static final String SUBSECTION_POSTERIORS
            = "<h3>Posterior probabilities for ancestral reconstruction</h3>" +
            "<p>Ancestral reconstruction " +
            "can be carried out by " +
            "computing posterior probabilities for the family sizes at " +
            "inner nodes. Let &xi;[<var>u</var>] denote the family size at " +
            "node <var>u</var>. The vector " +
            "&xi;=(&xi;[<var>u</var>]: all nodes <var>u</var>) " +
            "is a so-called phylogenetic character. The " +
            "phylogenetic birth-and-death model defines the distribution " +
            "of &xi;. For a phylogenetic profile &Phi;, " +
            "Count computes posterior probabilities for different " +
            "historical characteristics conditioned on " +
            "{&xi;[<var>L</var>]=&Phi;} " +
            "= {&xi[<var>u</var>]=&Phi;[<var>u</var>] for all leaves <var>u</var>&isin;<var>L</var>}." +
            "In particular, the following probabilities are computed.</p>" +
            "<ul>" +
            "<li>presence of at least one family member at an ancestral node <var>u</var></li>: " +
            "   p<sub>&Phi;</sub>(<var>u</var>:&ge;1) = Pr{&xi;[<var>u</var>]&ge;1 | &xi;[<var>L</var>]=&Phi;}" +
            "<li>multiple members at an ancestral node <var>u</var>: " +
            "   p<sub>&Phi;</sub>(<var>u</var>:<tt>m</tt>)=Pr{&xi;[<var>u</var>]&ge;2 | &xi;[<var>L</var>]=&Phi;}</li>" +
            "<li>family gain (size increase from 0 to a positive value) on an edge <var>u</var><var>v</var>: " +
            "   p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>g</tt>) = Pr{&xi;[<var>u</var>]=0, &xi;[<var>v</var>]&ge;1 | &xi;[<var>L</var>]=&Phi;}</li>" +
            "<li>family loss (size decrease from a positive value to 0) on an edge <var>u</var><var>v</var>: " +
            "   p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>l</tt>) = Pr{&xi;[<var>u</var>]&ge;1, &xi;[<var>v</var>]=0 | &xi;[<var>L</var>]=&Phi;}</li>" +
            "<li>family expansion (size increase from 1 to a larger value) on an edge <var>u</var><var>v</var>: " +
            "   p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>++</tt>) = Pr{&xi;[<var>u</var>]=1, &xi;[<var>v</var>]&ge;2 | &xi;[<var>L</var>]=&Phi;}</li>" +
            "<li>family contraction (size decrease from a larger value to 1) on an edge <var>u</var><var>v</var>: " +
            "   p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>--</tt>) = Pr{&xi;[<var>u</var>]&ge;2, &xi;[<var>v</var>]=1 | &xi;[<var>L</var>]=&Phi;}</li>" +
            "</ul>" +
            "<h3>Expectations</h3>"+
            "<p>A great advantage of using " +
            "posterior probabilities is that they can be summed together " +
            "to obtain expectations, which are excellent aggregate " +
            "characteristics for family dynamics and ancestral lineages.</p>" +
            "<p>Count uses five characteristics for family dynamics, " +
            "computed by summing posteriors probabilities along all edges for the same family.</p>" +
            "<ul>" +
            "<li>number of family gains across all lineages for a family <var>f</var> with profile &Phi;<sub><var>f</var></sub>=&Phi;: " +
            "   E<sub><var>f</var></sub>(<tt>g</tt>) = &sum;<sub><var>u</var><var>v</var></sub> " +
            "       p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>g</tt>)</li>" +
            "<li>number of family losses across all lineages for a family <var>f</var> with profile &Phi;<sub><var>f</var></sub>=&Phi;: " +
            "   E<sub><var>f</var></sub>(<tt>l</tt>) = &sum;<sub><var>u</var><var>v</var></sub> " +
            "       p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>l</tt>)</li>" +
            "<li>number of family expansions across all lineages for a family <var>f</var> with profile &Phi;<sub><var>f</var></sub>=&Phi;: " +
            "   E<sub><var>f</var></sub>(<tt>++</tt>) = &sum;<sub><var>u</var><var>v</var></sub> " +
            "       p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>++</tt>)</li>" +
            "<li>number of family contractions  across all lineages for a family <var>f</var> with profile &Phi;<sub><var>f</var></sub>=&Phi;: " +
            "   E<sub><var>f</var></sub>(<tt>g</tt>) = &sum;<sub><var>u</var><var>v</var></sub> " +
            "       p<sub>&Phi;</sub>(<var>u</var><var>v</var>:<tt>--</tt>)</li>" +
            "<li>number of arrivals: presence at root + gains on all edges for a family <var>f</var> with profile &Phi;<sub><var>f</var></sub>=&Phi;: " +
            "   E<sub><var>f</var></sub>(<tt>A</tt>) = p<sub>&Phi;</sub>(root:&ge;1)" +
            "       + E<sub><var>f</var></sub>(<tt>g</tt>)</li>"+
            "</ul>" +
            "<p>By summing across all families for an ancestral node or an edge, " +
            "Count computes expectations to quantify changes in different lineages." +
            "<ul>" +
            "<li>number of families present at an ancestral node <var>u</var>: " +
            "   E(<var>u</var>:&ge;1) = &sum;<sub><var>f</var></sub> p<sub>&Phi;<sub><var>f</var></sub></sub>(<var>u</var>:&ge;1) </li>" +
            "<li>number of multi-member families present at an ancestral node <var>u</var>: " +
            "   E(<var>u</var>:<tt>m</tt>) = &sum;<sub><var>f</var></sub> " +
            "       p<sub>&Phi;<sub><var>f</var></sub></sub>(<var>u</var>:<tt>m</tt>) </li>" +
            "<li>number of families gained on an edge <var>u</var><var>v</var>: " +
            "   E(<var>u</var><var>v</var>:<tt>g</tt>) = &sum;<sub><var>f</var></sub> " +
            "       p<sub>&Phi;<sub><var>f</var></sub></sub>(<var>u</var><var>v</var>:<tt>g</tt>) </li>" +
            "<li>number of families lost on an edge <var>u</var><var>v</var>: " +
            "   E(<var>u</var><var>v</var>:<tt>l</tt>) = &sum;<sub><var>f</var></sub> " +
            "       p<sub>&Phi;<sub><var>f</var></sub></sub>(<var>u</var><var>v</var>:<tt>l</tt>) </li>" +
            "<li>number of families expanded on an edge <var>u</var><var>v</var>: " +
            "   E(<var>u</var><var>v</var>:<tt>++</tt>) = &sum;<sub><var>f</var></sub> " +
            "       p<sub>&Phi;<sub><var>f</var></sub></sub>(<var>u</var><var>v</var>:<tt>++</tt>) </li>" +
            "<li>number of families contracted on an edge <var>u</var><var>v</var>: " +
            "   E(<var>u</var><var>v</var>:<tt>--</tt>) = &sum;<sub><var>f</var></sub> " +
            "       p<sub>&Phi;<sub><var>f</var></sub></sub>(<var>u</var><var>v</var>:<tt>--</tt>) </li>" +
            "</ul>" +
            "";
    private static final String SUBSECTION_ABSENT_FAMILIES_T = "Absent families";
    private static final String SUBSECTION_ABSENT_FAMILIES
            = "<p>An <em>absent family</em> is a family that has a phylogenetic " +
            "profile &Phi;=(0,0,...,0), i.e., one that has no " +
            "members at any terminal node. Absent families are immaterial " +
            "in parsimony analyses, but the phylogenetic birth-and-death " +
            "model assigns a well-defined probability <var>p</var><sub>0</sub> " +
            "to the all-0 profile. The likelihood optimization assumes that the data " +
            "set does not include any families with an all-0 profile, and " +
            "corrects the likelihood formula appropriately. Even for an " +
            "absent family, there is a small probability that the family " +
            "history includes at least one ancestral presence. " +
            "If <var>p</var><sub>0</sub> is the probability of an all-0 profile, " +
            "then the number of absent families is estimated as " +
            "<var>n</var><sub>0</sub> = <var>n</var><var>p</var><sub>0</sub>/(1-<var>p</var><sub>0</sub>). " +
            "Count can take absent families into account when " +
            "inferring lineage-specific statistics. " +
            "For instance, the number of family gains on an edge is corrected by the " +
            "number of gains in all-0 profiles: " +
            "E'(<var>u</var><var>v</var>:<tt>g</tt>) = E(<var>u</var><var>v</var>:<tt>g</tt>) " +
            "   + <var>n</var><sub>0</sub> p<sub>(0,...,0)</sub>(<var>u</var><var>v</var>:<tt>g</tt>). " +
            "Similar corrections are applied to every statistic computed across all families. "+
            "</p>";

    private static final String CHAPTER_USING_COUNT_T = "Using Count";
    private static final String SECTION_DESIGN_CONCEPTS_T = "Design concepts";
    private static final String SUBSECTION_ERRORS_T = "Errors";
    private static final String SUBSECTION_ERRORS
            = "<p>It may happen that something goes wrong. Count displays a " +
            "window with the error message in such cases. If you think " +
            "that the error was caused by a programming bug, or you would " +
            "like to ask me for help with the situation, then send me an email that includes " +
            "the detailed technical error message in the body. " +
            "The technical details are shown only when you select " +
            "the corresponding checkbox.</p>";
    public static final String SUBSECTION_SESSIONS_WORKAREA_T = "Sessions and the work area";
    public static final String SUBSECTION_SESSIONS_WORKAREA
            = "<p>Count operates with <em>sessions</em>: each session is associated with a " +
            "fixed species phylogeny. More than one session may be open at one time: e.g., " +
            "the same data set may be analyzed with different phylogenies simultaneously. " +
            "A session has three main components, represented by the tabs of the displayed " +
            "workspace: a species phylogeny (Tree), " +
            "a <em>browser</em> for data sets and analysis results " +
            "(Data) and a browser for probabilistic models (Rates)." +
            "The current session's name is displayed as the window title.</p>";

    public static final String SUBSECTION_BROWSERS_T = "Browsers, primary items and views, exporting";
    public static final String SUBSECTION_BROWSERS
            = "<p>The Data and Rates tabs are attached to " +
            "browser displays. A browser consists of a hierarchy on the left," +
            "and an information panel on the right, corresponding to the item selected " +
            "in the hierarchy. " +
            "<em>Primary items</em> in the hierarchy (depicted as folders of a file system)" +
            "are alignments or data tables (under the Data tab), or " +
            "rate models (under the Rates tab). " +
            "Primary items may have <em>views</em>, which correspond to various " +
            "analysis tasks. Views are descendant nodes in the hierarchy (depicted as " +
            "documents or bullets, depending on the operating system).</p>" +
            "<p>Nodes of the hierarchy have small " +
            "associated popup menus which you can bring up by right-clicking on them " +
            "(or by Ctrl-click on a Mac). The popup menu items include the removal of the node  " +
            "from the browser, and possibly saving options. Intron tables and rate models can be saved, " +
            "and views are typically <em>exported</em> " +
            "into text files. (The difference is that exported views cannot be loaded later, " +
            "but saved tables and rate models can.)</p>" +
            "<p>The browsers operate with split panes. You can expand a pane completely by clicking on " +
            "the little triangles on the bottom or the far right of the dividers, or resize the panes " +
            "by dragging the dividers with the mouse. </p>";

    public static final String SUBSECTION_TREE_DISPLAYS_T = "Tree displays";
    public static final String SUBSECTION_TREE_DISPLAYS
            = "<p>There are several graphical displays that show analysis results on the " +
            "session's phylogeny. Tree displays have some associated control elements " +
            "in the bottom tool bar. Most importantly, there is a zooming spinner on the " +
            "bottom right, where you can set a relative magnification " +
            "factor (displayed as a percentage). " +
            "You can select a tree node by clicking on it. In that way, you " +
            "get some more specific information about the selected node: " +
            "it depends on the context what that information exactly is. " +
            "You can select at most one node at a time. " +
            "In order to deselect all nodes, click somewhere away " +
            "from the tree nodes within the tree display.</p>";
    public static final String SUBSECTION_TABLE_COLUMNS_T = "Table displays: column rearrangements and row sorting";
    public static final String SUBSECTION_TABLE_COLUMNS
            = "<P>Most results are shown in tree displays together " +
            "with table displays. Row selection in the tables affects the " +
            "tree display and node selections in the tree display may " +
            "affect the row selection in a table for lineages.</p>" +
            "<p>Columns " +
            "can be rearranged at will by dragging the column headers.</p>" +
            "<p>Table displays can be sorted row-wise by clicking on a column " +
            "header. The column by which the table is sorted has a mark " +
            "next to it, also indicating the direction (increasing or " +
            "decreasing order).</p>" +
            "<p>Count " +
            "works with high-precision numerical values internally " +
            "(Java's 64-bit floating-point type <tt>double</tt>), " +
            "but the tables use rounding. Zero is " +
            "denoted by a dot. The cell's tool tip gives the exact value." +
            "</p>";

    public static final String SUBSECTION_COMMENTS_T = "Comments in input and output files";
    public static final String SUBSECTION_COMMENTS
            ="<p>Output and input files may contain comments in lines starting with " +
            "<tt>#</tt>. Such lines are ignored on input, and do not " +
            "contain essential information. If necessary (e.g., in order " +
            "to prepare for import into Excel), they are easily filtered " +
            "out in Unix on the command-line.</p>" +
            "\n<pre>\n" +
            "% grep -v '#' file > stripped\n" +
            "</pre>\n";
   public static final String SECTION_SESSIONS_T = "Sessions";
   public static final String SUBSECTION_PHYLOGENY_T = "Phylogeny";
   public static final String SUBSECTION_PHYLOGENY 
           = "<p>Data analysis in Count starts with opening a new session " +
           " (Menu: <b>Session-&gt;Start new session...</b>). " +
           "A session is opened by loading a species phylogeny. The phylogeny is expected " +
           "to be in Newick format (<a href=\"http://evolution.genetics.washington.edu/phylip/newicktree.html\">http://evolution.genetics.washington.edu/phylip/newicktree.html</a>)" +
           "used by Phylip and other fine software packages for molecular evolution. " +
           "The branch lengths of this phylogeny are ignored in all cases, except " +
           "when computing Propensity for Gene Loss (PGL). The inner nodes of " +
           "the tree may have more than two children: Count can deal with " +
           "arbitrary multifurcations.</p>" +
           "<p>The phylogeny is displayed under the Tree tab. For " +
           "convenience in the graphical user interface, it is " +
           "recommended that you use short names (3-4 letters) for the " +
           "terminal taxa. The inner nodes of the tree are numbered as " +
           "1,2,... (in a postorder traversal). It is very useful to " +
           "name the inner nodes of the phylogeny, which is possible in " +
           "Newick format:</p>\n"+
           "<pre>\n"+
           "((Natph, Halsp, Halwa) Halobacteriales,\n" +
           " ((Metcu, Methu, Metla) Methanomicrobiales,\n"+
           "  (Metbu,Metsa,(Metac, Metma, Metba) Methanosarcina) Methanosarcinales \n"+
           " ) Methanomicrobia \n"+
           ") root;\n" +
           "</pre>"+
           "<p><strong>Attention:</strong> in the Newick format specification, underscore characters " +
           "appearing in taxon names are to be converted into whitespace. For instance, <tt>E_coli</tt> " +
           "in the tree file will be interpreted as \"E coli\". This is important to take into account for naming " +
           "the corresponding columns in the family size table, because the converted name is used there. " +
           " In the example, the table should have a column named \"E coli\" and not \"E_coli.\"</p>"+
           "<p>The open sessions are listed under the <b>Session</b> menu, " +
           "where they can be selected to switch back and forth between " +
           "them. There is a pointing hand icon next to the current active section " +
           "(shown in the window title). The active session can be " +
           "disposed of by closing it (Menu: <b>Session&gt;Close session</b>). " +
           "</p>";
   public static final String SUBSECTION_SAVING_T = "Saving your work";
   public static final String SUBSECTION_SAVING
           ="<p>It is possible to save all your sessions, together " +
           "with all the rate models, data tables, and analysis views at once: " +
           "<b>Session&gt;Save everything ...</b>. " +
           "The saved file is in a machine-readable format (XML), " +
           "and can be loaded later to restore your analysis pipeline:  " +
           "<b>Session&gt;Open previously saved session(s) ...</b>. " +
           "Note that the menu point is only available when there are no " +
           "sessions open yet. The saved XML file can be compressed with Gzip " +
           "(<tt>gzip saved-sessions.xml</tt> in the command line), " +
           "and loaded later directly: Count checks for a <tt>.gz</tt> " +
           "extension, and uncompresses the file on the fly.</p>";
   public static final String SUBSECTION_TABLE_T = "Family size table";
   public static final String SUBSECTION_TABLE
           = "<p>The input data, on which various analyses can be performed, " +
           "is a family size <strong>table</strong>. " +
           "The table of phylogenetic profiles is a TAB-delimited text " +
           "file. Every row corresponds to a homolog gene family, with " +
           "the exception of the first row that gives the column headers. " +
           "The first column is the family name. The second, third, etc. " +
           "columns correspond to terminal taxa of the phylogenetic tree. The " +
           "column headers must specify the terminal taxon names, in an " +
           "arbitrary order. Columns with taxons missing from the " +
           "phylogeny are ignored. </p>\n" +
           "<pre>\n" +
           "family	Aerpe	Arcfu	Calma	Censy	Halma	Halsp\n"+
           "arCOG00001	1	2	1	0	0	0\n" +
           "arCOG00002	1	0	1	1	1	1\n" +
           "arCOG00004	0	0	0	0	1	1\n"+
           "</pre>\n" +
           "<p>Note that the relevant files of the COG, KOG and arCOG " +
           "databases can be used immediately, without any input format conversion.</p>"+
           "<p>The parsimony analyses can handle missing entries in the "+
           "table, denoted by <tt>?</tt>, but other programs treat " +
           "missing data as a family size of 0.</p>" +
           "<p>A family size table can be opened from the menu: " +
           "<b>Data&gt;Open table ...</b>. The opened table is displayed " +
           "in the Data tab's browser as a primary item. The displayed " +
           "columns include <tt>#lin</tt> and <tt>#mem</tt>, " +
           "which are the total number of terminal taxa " +
           "with at least one member, and total number of family " +
           "members, respectively. Other columns are the family indices " +
           "(original order in the file), family names, and number of " +
           "homologs per terminal taxon (column <tt>&Phi;</tt><var>x</var> for taxon " +
           "<var>x</var>). </p>";
   public static final String SUBSECTION_ANNOTATIONS_T = "Family annotations";
   public static final String SUBSECTION_ANNOTATIONS
           = "<p><em>Annotations</em> " +
           "define family properties: the default annotation is simply the family " +
           "name. Further annotations (e.g., COG functional category) " +
           "can also be included in Count. Note that annotations are " +
           "text fields: numerical annotations are not supported. There " +
           "are two ways to annotate families:  an annotated family size table can be " +
           "opened (<b>Data&gt;Open annotated table...</b>), or family annotations can be " +
           "loaded from a separate file (<b>Data&gt;Load family annotations...</b>).</p>"+
           "<p>If you load an annotated table (<b>Data&gt;Open annotated table...</b>), " +
           "then every column with a header not correponding to a terminal taxon name " +
           "is considered as an annotation column. You can have, for instance, <q>Category</q> " +
           "and <q>Description</q> columns. Column headers must be unique " +
           "(you cannot have two <q>Category</q> columns). You can also " +
           "load a simple family size table. In that particular case, " +
           "there will be an annotation column for every organism that " +
           "is not present at the terminal taxa of the session's " +
           "phylogeny, although family sizes there will treated as text " +
           "values. Family name is always taken from the first column of " +
           "the input file.</p>"+
           "<p>You can also add annotations from a separate file at any " +
           "time (<b>Data&gt;Load family annotations...</b>). " +
           "The annotations file is a text table, where the fields " +
           "can be separated by comma (<tt>.csv</tt> extension), TAB " +
           "(usually <tt>.txt</tt> extension), " +
           "or any other character. After the file is selected, the text format, " +
           "and the imported columns are selected through a popup dialog. " +
           "Note that COG family definition files (such as <tt>arCOGdef.csv</tt>) can be " +
           "used immediately. " +
           "The first column of the annotation file must give the family " +
           "name: families can be listed in arbitrary order. You can " +
           "select the annotation columns that you would like to import " +
           "through the dialog. The displayed column headers can be " +
           "edited. They are either SKIPPED, meaning that they will not " +
           "be imported, or have a different title from <q>SKIPPED.</q> " +
           "There are two predefined options: <q>Category</q> and " +
           "<q>Description,</q> but you can select any other text, as long " +
           "as the columns have different names. The imported annotation " +
           "columns are added to the currently selected table display, " +
           "and all its descendant views.</p>";
   public final static String SUBSECTION_FILTERING_T = "Family selections, filtering, and " +
           "absence/presence transformations";
   public final static String SUBSECTION_FILTERING 
           = "<p>In display tables for family size and analysis results, you " +
           "can select multiple families using the mouse directly, or by " +
           "using logical selection criteria. Selection criteria are " +
           "displayed by double-clicking on a table cell, in a popup " +
           "menu. If you double-click on a numerical column, then the " +
           "selection options are \"equal,\" \"less than or equal to,\" " +
           "and \"greater than or equal to,\" with the reference value " +
           "taken from the cell you clicked on. If you double-click on a " +
           "text column (annotations and family name), then the " +
           "selection options are \"equals\" and \"contains.\" In this " +
           "way, you can select families with a particular functional " +
           "category, or size, or taxon representation. Or, using the " +
           "table displaying some analysis results, you can define " +
           "selection criteria based on presence at ancestral nodes, or " +
           "other inferred characteristics.</p>" +
           "<p>Selected families can be extracted into a new table " +
           "(<b>Data&gt;Extract selected families into a new table</b>). " +
           "The filtered table appears as a descendant view in the " +
           "Data browser.<p>" +
           "<p>Finally, family sizes can be transformed into binary " +
           "profiles (<b>Data&gt;Transform numerical profiles into binary (presence/absence) profiles</b>). " +
           "Every positive value is replaced by '1' in the result. The binary " +
           "table is displayed as a descendant view in the Data browser.</p>" +
           "<p> You can save the tables you created through filtering and binary " +
           "transformations. Family size tables can be saved through the " +
           "popup menu for the corresponding node in the data browser.</p>";
   public static final String SECTION_RATES_T = "Rates";
   private static final String SECTION_RATES
           = "<p>Count uses phylogenetic birth-and-death " +
           "models in probabilistic inference. Model parameters (rates) can be set by optimization " +
           "on the currently selected family size table, or previously computed rate models can be " +
           "loaded from a <em>rate model file</em>. Rate models are" +
           "displayed as <em>rate panels</em> under the Rates tab.</p>";
   public static final String SUBSECTION_RATE_FILE_T = "Rate model file";
   private static final String SUBSECTION_RATE_FILE
           = "<p>Model parameters are given in a" +
           "text file. The simplest way of " +
           "browsing the rate file is to strip comments " +
           "(<tt>grep -v '#'</tt>) which gives a TAB-delimited table " +
           "that can be imported into Excel or other spreadsheet program. " +
           "The table columns " +
           "are the edge length <var>t</var><sub><var>e</var></sub>, " +
           "duplication rate &lambda;<sub><var>e</var></sub>, " +
           "loss rate &mu;<sub><var>e</var></sub>, and " +
           "gain rate &kappa;<sub><var>e</var></sub>, " +
           "followed by columns with debugging information. Every row " +
           "corresponds to an edge in the phylogeny, enumerated in a " +
           "postorder traversal. (The debug information starts with the " +
           "name of the edge's child node.) The rates and edge lengths " +
           "are scaled by the total gene loss rate so that " +
           "&mu;<sub><var>e</var></sub>=1 on every edge <var>e</var>.</p>\n" +
           "<pre>\n" +
           "length (t)    duplication (lambda)    loss (mu)   transfer (kappa)  // ...\n"+
           "1.5820017453696325    0.38702704740849914 1.0  0.012376046605348303  // ...\n" +
           "4.139742169467891    0.08598802637984121  1.0  0.0021845554103576185 // ...\n" +
           "0.23206654225569434   0.47774576195966995  1.0 0.016986771452090894  // ... \n" +
           "0.157023248086316    0.569613362320818    1.0  0.012144222758613127  // ...\n" +
           "</pre>\n" +
           "<p>The last lines of the rates file give the remaining model parameters such as " +
           "the distribution of family-specific rate factors " +
           "<var>t</var>'<sub><var>f</var></sub>, " +
           "&lambda;'<sub><var>f</var></sub>, "+
           "&mu;'<sub><var>f</var></sub>, "+
           "&kappa;'<sub><var>f</var></sub>, and the family size " +
           "distribution at the root. In the example below, gene loss " +
           "&mu;'<sub><var>f</var></sub> (\"loss\") " +
           "and gain &kappa;<sub><var>f</var></sub> (\"transfer\") " +
           "are constant, but " +
           "gene duplication &lambda;'<sub><var>f</var></sub> (\"duplication\") " +
           "and duration <var>t</var>'<sub><var>f</var></sub> (\"length\") " +
           "have Gamma distributions (shape paremeters of 0.846... and 0.828...) " +
           "discretized using 4 categories. The family size at the root " +
           "has a Poisson distribution with mean 0.022...</p>\n" +
           "<pre>\n" +
           "|variation      duplication     4      0.8461581445706693      0.0 \n" +
           "|variation      loss    1      1.0     0.0 \n" +
           "|variation      transfer        1       1.0    0.0 \n" +
           "|variation      length  4       0.8281017740809735 \n" +
           "|root  Poisson 0.02253828154584124\n" +
           "</pre>";
   public static final String SUBSECTION_RATES_PANEL_T = "Rates panel";
   private static final String SUBSECTION_RATES_PANEL
           = "<p>The information panel for a rate " +
           "model consists of three parts: a table showing numerical " +
           "values of gain/loss/duplication rates (on the upper left), a " +
           "graphical illustration of rate categories (on the upper " +
           "right), and a graphical illustration of branch-specific " +
           "gain, loss and duplication rates (in the lower half).</p>" +
           "<h3>Rates panel: table</h3>" +
           "<p>The table on the upper left-hand side gives the branch-specific " +
           "prototypical gain, loss and duplication rates " +
           "(&kappa;<sub><var>e</var></sub> <var>t</var><sub><var>e</var></sub>, " +
           "&mu;<sub><var>e</var></sub> <var>t</var><sub><var>e</var></sub>, " +
           "and &lambda;<sub><var>e</var></sub> <var>t</var><sub><var>e</var></sub>, " +
           "respectively) for each branch <var>e</var>. " +
           "Branches are specified by the nodes they lead to. There are " +
           "no rates next to the root node. If a table row is selected, " +
           "then the corresponding tree node " +
           "is highlighted in the tree display on the bottom.</p>" +
           "<h3>Rates panel: graphical display of rate variation</h3>" +
           "<p>The prior distributions of family-specific rate factors " +
           "are defined by possible no-gain (&kappa;'<sub><var>f</var></sub>=0) " +
           "and no-duplication categories (&lambda;'<sub><var>f</var></sub>=0), " +
           "and possible rate factors for the discretized Gamma distributions " +
           "in case of edge length (<var>t</var>'<sub><var>f</var></sub>), " +
           "loss rate (&mu;'<sub><var>f</var></sub>), " +
           "duplication rate (&lambda;'<sub><var>f</var></sub>) " +
           "and gain rate (&kappa;'<sub><var>f</var></sub>). " +
           "The Gamma distribution plots also give the shape " +
           "parameter &alpha;, and shade the corresponding continuous distribution.</p>" +
           "<h3>Rates panel: tree display</h3>" +
           "<p>The bottom part of the rates panel displays branch-specific " +
           "model components (&kappa;<sub><var>e</var></sub> <var>t</var><sub><var>e</var></sub>, " +
           "&mu;<sub><var>e</var></sub> <var>t</var><sub><var>e</var></sub>, " +
           "and &lambda;<sub><var>e</var></sub> <var>t</var><sub><var>e</var></sub>). " +
           "The tree panel also shows the prior family size distribution at the root. " +
           "Rates may vary much across different lineages, and therefore " +
           "it is not possible to have proportional branch lengths everywhere in " +
           "the display. Instead, an \"informative ellipsis\" is used " +
           "for long branches: the ratio between the solid part of the " +
           "branch and the entire plotted branch length equals the ratio " +
           "of the displayed branch length and the true branch length. " +
           "For instance, if the displayed branch length corresponds to " +
           "a loss rate of &lambda;=&lambda;<sub><var>e</var></sub> <var>t</var><sub><var>e</var></sub>=1, " +
           "and the true loss rate is &lambda=4;, " +
           "then one-quarter of the displayed branch is solid and the rest is dotted.</p>" +
           "<p>A legend panel is laid over the tree panel on the left. The legend " +
           "panel can be disabled by clicking the <b>Legend</b> checkbox in the " +
           "bottom tool bar. Other checkboxes (<b>Loss</b>, <b>Duplication</b>, <b>Gain</b>) " +
           "are used to select the rate components that are shown in the " +
           "tree panel and the legend. The legend panel shows the " +
           "scaling for the different rate components at the actual zoom " +
           "level (set by the the bottom tool bar's spinner on the " +
           "right).</p>" +
           "<p>Additional information is shown for the selected tree node " +
           "(selection is done either by clicking on it directly in the " +
           "tree panel, or by clicking on the corresponding row of the " +
           "table on the upper left). In particular, the distributions " +
           "for inparalog and xenolog group sizes are shown. The " +
           "distribution plots show the probabilities for group sizes 0, " +
           "1, 2, ...; the bars are scaled linearly so that the Y " +
           "axis is of length 1.</p>";
   public static final String SUBSECTION_RATE_OPTIMIZATION_T = "Rate model optimization";
   private static final String SUBSECTION_RATE_OPTIMIZATION
           = "<p>Count computes the " +
           "model parameters of the phylogenetic birth-and-death model " +
           "by the numerical optimization of the likelihood. " +
           "The likelihood computation assumes that " +
           "there are no all-0 profiles in the " +
           "data set. It is therefore recommended that you first filter " +
           "those families out before optimizing the likelihood. The " +
           "simplest way to do that is to sort the family table by " +
           "lineage-weight of the profile (column <tt>#lin</tt>), " +
           "and double-click on a cell with <tt>#lin</tt>=1. The popup " +
           "selection includes the option of <tt>#lin</tt>&ge; 1. " +
           "Selected families can be then filtered into a separate table " +
           "(<b>Data&gt;Extract selected families...</b>.</p>" +
           "<p>Prior to proceeding to the actual computation, optimization " +
           "parameters need to be set in the window that appears after " +
           "selecting the menu point <b>Rates&gt;Optimize rates...</b>. " +
           "Optimization parameters are grouped under the tabs " +
           "Model type and Model parameters.</p>" +
           "<h3>Rate optimizaton: model type</h3>" +
           "<p>First, the initial model needs to be selected: this can be Count's " +
           "predefined null model, or a previously computed rate model. " +
           "This latter option is offered only if a rate model is " +
           "selected in the Rates panel. The optimized model and " +
           "parameters are initialized using the selected model.</p>" +
           "<p>Second, the optimized model architecture needs to be " +
           "selected: gain-loss-duplication, duplication-loss, " +
           "gain-loss, and pure loss. The most general model is the " +
           "gain-loss-duplication model, where there is no restriction " +
           "on the lineage-specific rates. In a duplication-loss model, " +
           "all gain rates are zero (&kappa;<sub><var>e</var></sub>=0); " +
           "in a gain-loss model, all duplication rates are zero (&lambda;<sub><var>e</var></sub>=0). " +
           "In a pure loss model, both gain and duplication rates are zero.</p>" +
           "<p>Third, the type of the prior distribution at the root needs " +
           "to be selected: this may be Poisson, negative binomial, or " +
           "Bernoulli (point) distribution.</p>" +
           "<p>Fourth, it must be selected if duplication and gain rates " +
           "may differ between tree edges. If, say, the <b>Same gain/loss ratio in all lineages</b> " +
           "checkbox is selected, then the optimization assumes " +
           "that &kappa;<sub><var>e</var></sub>=&kappa; " +
           "for some common gain rate &kappa;, and optimizes the single " +
           "model parameter &kappa;.</p>" +
           "<p>Fifth, the type of the rate variation across families needs " +
           "to be chosen: this includes the number of discrete Gamma " +
           "categories (=1 if there is no Gamma variation), and " +
           "possible no-duplication and no-gain categories.</p>" +
           "<p>The final set of parameters comprises computational " +
           "parameters for the numerical optimization. The optimization " +
           "proceeds in rounds: all model parameters are optimized once " +
           "in each round. The optimization stops after the given " +
           "maximum of optimization rounds, or earlier, when in two " +
           "consecutive rounds, the log-likelihood (natural logarithm) " +
           "changes by less than the given convergence threshold.</p>" +
           "<h3>Rate optimization: model parameters</h3>" +
           "<p>Under the Model parameters tab, you can set the initial " +
           "values for all model parameters, as well as exclude certain " +
           "parameters from the optimization. In order to exclude some " +
           "parameter from the optimization, select its " +
           "<b>Fixed</b> checkbox.</p>" +
           "<p>The model parameters include the prior family size " +
           "distribution at the root (one or two parameters for Poisson, " +
           "negative binomial, or Bernoulli distribution), parameters " +
           "for rate variation across families (the set of tunable " +
           "parameters depends on the rate variation type selected under " +
           "the Model type tab), and lineage-specific " +
           "rates. Lineage-specific rates can be fixed individually, or all at " +
           "once using the \"master\" checkboxes in the \"all edges\" row.</p>" +
           "<p>For technical reasons, the rates and edge lengths are " +
           "scaled in such a way that the loss rate equals 1 on every " +
           "edge. Note that the edge lengths of the input phylogeny are " +
           "ignored in the probabilistic inference.</p>" +
           "<h3>Rate optimization: computing</h3>" +
           "<p>The actual optimization process starts when the " +
           "<b>Perform optimization</b> button is pressed. The progress " +
           "of the optimization can be followed in the rate model " +
           "display that shows up.</p>" +
           "<p>The optimization is launched in a background process, and " +
           "you can continue working with Count, performing other " +
           "analysis steps. The rate model display will be updated " +
           "continuously in the course of the optimization process. The " +
           "display for the optimized model appears in the Rates " +
           "browser. The bottom tool bar in optimized rate model " +
           "displays includes progress indicators: a progress bar " +
           "showing the current round, information about the current " +
           "optimization step, the value of the log-likelihood " +
           "(LL), and its increase in the last round " +
           "(&Delta;LL). You cannot use the rate model during " +
           "optimization, and you should not save it (because it changes " +
           "during the save). Instead, you can take <em>snapshots</em> " +
           "of the current rate model " +
           "during optimization, which will appear as descendant nodes " +
           "of the optimized model display in the Rates browser. " +
           "You can save the snapshot into a file, or perform ancestral " +
           "reconstruction with it. The bottom tool bar has two specific " +
           "buttons: one for taking snapshots of the current rate model " +
           "during optimization, and another button for canceling the " +
           "process (<b>Stop</b>).</p>" +
           "<p>After the optimization finished, you can use the optimized " +
           "rate model as any other rate model: you can save it, or " +
           "perform analyses with it. " +
           "Rate models (including snapshots and finished optimizations) can be " +
           "saved through the popup menu for the corresponding node in " +
           "the Rates browser.</p>" +
           "";
   public static final String SECTION_ANALYSIS_T = "Analysis tasks";
   private static final String SECTION_ANALYSIS
           = "<p>You can perform ancestral " +
           "inference and analyze family dynamics using the options " +
           "available under the <b>Analysis</b> menu point. " +
           "Namely, the following methods are implemented.</p>" +
           "<ul>" +
           "<li>Dollo parsimony</li>" +
           "<li>Wagner parsimony</li>" +
           "<li>posterior probabilities by the selected phylogenetic birth-and-death model</li>" +
           "<li>Propensity for Gene Loss (PGL)</li>" +
           "</ul>";
   public static final String SUBSECTION_ANALYSIS_PANEL_T = "Analysis panels";
   private static final String SUBSECTION_ANALYSIS_PANEL
           = "<p>The analysis panels have the same basic design " +
           "concept: they all consist of three parts. The three parts " +
           "are (1) a table for family-specific information on the upper " +
           "left, (2) a table for lineage-specific information on the " +
           "upper right, and (3) a tree display on the bottom.</p>" +
           "<h3>Analysis panel: family table</h3>" +
           "<p>The family table on the upper left has a row for each " +
           "family. Its columns include family index, family name, " +
           "possibly family annotations, number of terminal lineages " +
           "the family is present in (\"#lin\"), total number of " +
           "members at terminal lineages (\"#mem\"), and " +
           "phylogenetic profile. Additional columns are specific to the " +
           "analysis methods.</p>" +
           "<p>The profile is depicted graphically: " +
           "black bars show presence, " +
           "with a height that is proportional to the logarithm of the " +
           "family size at each node. The tool tip for a profile cell " +
           "gives the exact numerical profile. </p>" +
           "<p>Multiple rows can be selected in the family table in the " +
           "usual manner for your operating system (e.g., shift+click " +
           "for range selection, Cmd-A or Ctrl-A for selecting all rows, " +
           "etc.). The lineage table on the upper right shows sums (of " +
           "gains, losses etc.) across the selected families by " +
           "lineages. The tree display on the bottom illustrates the " +
           "history of the selected families.</p>" +
           "<h3>Analysis panel: lineage table</h3>" +
           "<p>The lineage table " +
           "on the upper right gives aggregate information over the " +
           "selected families. Table rows correspond to lineages.</p>" +
           "<p>Table columns may include total number of families present " +
           "(\"Families\") and total number of multi-member " +
           "families (<tt>:m</tt>) present at the node, as well as event " +
           "totals on the edge leading to the node: family gains " +
           "(<tt>:g</tt>), family losses (<tt>:l</tt>), expansions " +
           "(<tt>++</tt>), and contractions (<tt>--</tt>).</p>" +
           "<h3>Analysis panel: tree display</h3>" +
           "<p>The tree display " +
           "on the bottom illustrates the inferred history of the " +
           "selected families. The bottom tool bar in analysis panels " +
           "gives information about the current selection in the " +
           "family table. For less than seven selected families, " +
           "their presence at nodes is illustrated individually, " +
           "otherwise only aggregated values " +
           "are shown by horizontal bars. Nodes in the tree panel " +
           "can be selected directly, or through the lineage table. " +
           "There is additional information displayed at the selected node.</p>" +
           "<p>Empty rectangles denote the absence of a family " +
           "from an ancestral node, and shaded rectangles denote " +
           "presence. The top half gives information about the " +
           "existence of multiple members. When it makes sense, " +
           "the horizontal extent of the shading in individual " +
           "rectangles is proportional to the likelihood of " +
           "the family presence, and that of having multiple members.</p>" +
           "<p>In addition to the inferred composition with " +
           "respect to selected families, family gains and losses " +
           "are indicated by green and orange bars on the edge " +
           "leading to each node. The solid part " +
           "for the two bars shows net gain or loss. " +
           "Numerical values, including those for family " +
           "expansions and contractions, are shown only for the selected node." +
           "</p>";
   public static final String SUBSECTION_ANALYSIS_DOLLO_T = "Analysis: Dollo parsimony";
   private static final String SUBSECTION_ANALYSIS_DOLLO
           = "<p>Ancestral presence in Dollo parsimony is inferred " +
           "by assuming that each family appeared only once, and that " +
           "the presence-absence pattern is " +
           "explained by lineage-specific losses. " +
           "The ancestral reconstruction by Dollo parsimony is accessed " +
           "through the menu <b>Analysis&gt;Family history by Dollo parsimony</b>. " +
           "</p>" +
           "<p>Multiple members are ignored in Dollo parsimony. Therefore, " +
           "the lineage table does not give expansion and contraction " +
           "counts, and the tree display illustrates family-level " +
           "characteristics only (presence, gain and loss). " +
           "If individual families are shown, then the first appearance of " +
           "each family is indicated by rectangles with bold green " +
           "frame. Rectangles with bold orange frame indicate " +
           "lineage-specific losses. The detailed information about the " +
           "selected node lists family gain and loss events leading to " +
           "the node.</p>" +
           "<p>Analysis results can be exported into a TAB-delimited text " +
           "file through the popup menu for the corresponding node in " +
           "the Data browser.</p>";
   public static final String SUBSECTION_ANALYSIS_WAGNER_T = "Analysis: Wagner parsimony";
   public static final String SUBSECTION_ANALYSIS_WAGNER
           = "<p>Wagner parsimony " +
           "penalizes the loss and gain of individual family members, " +
           "and infers the history with the minimum penalty. " +
           "The ancestral reconstruction by " +
           "Wagner parsimony is accessed through the menu " +
           "<b>Analysis&gt;Family history by Wagner parsimony</b>. " +
           "</p>" +
           "<p>The bottom tool bar includes a spinner for setting the gain " +
           "penalty (loss penalty is always =1) in asymmetric Wagner parsimony. " +
           "Recomputing the history with a new gain penalty value may " +
           "take some time: the process is launched in the background.</p>" +
           "<p>At the selected tree node, there is detailed information " +
           "about individual families: presence/absence " +
           "multiple members (<tt>m</tt>), and " +
           "family events on the edge leading to the node.</p>" +
           "<p>Analysis results can be exported into a TAB-delimited text " +
           "file through the popup menu for the corresponding item in " +
           "the Data browser.</p>";
   public static final String SUBSECTION_ANALYSIS_POSTERIORS_T = "Analysis: Posteriors";
   public static final String SUBSECTION_ANALYSIS_POSTERIORS
           ="<p>Ancestral reconstruction " +
           "by posteriors is " +
           "accessed through the menu " +
           "<b>Analysis&gt;Family history by posterior probabilities</b>.</p>" +
           "<p>The family table on the upper left shows the computed " +
           "ancestral reconstruction statistics and " +
           "family dynamics. Hover with the mouse over " +
           "a cell to see the explanation of the displayed value.</p>" +
           "<p>The bottom tool bar " +
           "includes a checkbox for including absent families in the " +
           "lineage-specific statistics. Computing the " +
           "history may take substantial time (minutes or even hours): " +
           "the process is launched in the background, " +
           "and you can continue working with Count.</p>" +
           "<p>At the selected tree node, there is detailed information about " +
           "lineage-specific statistics. " +
           "When individual families are shown at the selected node, " +
           "then the presence probability is followed by the statistics " +
           "with non-negligible values, including multi-member families " +
           "(<tt>m</tt>), gains, losses, " +
           "expansions (<tt>++</tt>) and contractions (<tt>--</tt>). " +
           "Probabilities different from 0 and 1 are given in " +
           "parentheses.</p>" +
           "<p>You can extract the posterior reconstruction into a " +
           "TAB-delimited text file through the popup menu of the " +
           "corresponding item in the Data browser. After " +
           "selecting the file, you can specify which columns the file " +
           "should include. The choices include posterior probabilities " +
           "of rate categories, various statistics for the ancestral " +
           "reconstruction and the family dynamics, " +
           "as well as the family annotations. The output file will have " +
           "an additional row for the absent profile if the correspnding " +
           "checkbox is selected. The ouput file format is discussed in " +
           "more detail in the PDF version of the manual.</p>"; 
   public static final String SUBSECTION_ANALYSIS_PGL_T = "Analysis: Propensity for gene loss";
   public static final String SUBSECTION_ANALYSIS_PGL 
           = "<p>Count can compute the so-called PGL (propensity for gene loss) " +
           "index, introduced by Krylov et al. [2003]. " +
           "PGL computation is is accessed through the menu " +
           "<b>Analysis&gt;PGL: propensity for gene loss (Krylov-Wolf-Rogozin-Koonin)</b>. " +
           "PGL is defined for a family " +
           "as (&sum;<sub><var>e</var></sub> loss<sub><var>e</var></sub>) length<sub><var>e</var></sub>/" +
           "(&sum;<sub><var>e</var></sub> length<sub><var>e</var></sub>), " +
           "where length<sub><var>e</var></sub> " +
           "is edge length (e.g., time between speciation events measured " +
           "in million years), and loss<sub><var>e</var></sub> " +
           "is an indicator for the optimal Dollo " +
           "parsimony reconstruction: loss<sub><var>e</var></sub>=1 " +
           "if the reconstruction posits a loss on edge <var>e</var>, " +
           "otherwise loss<sub><var>e</var></sub>=0. " +
           "The summation goes over the edges in " +
           "the subtree rooted at the first appearance of the family. " +
           "Count uses the edge lengths in the session's main " +
           "phylogeny.</p>" +
           "<p>Typically, one is interested in families that originate " +
           "at the same ancestor. In Count, you can select all such " +
           "families by performing Dollo parsimony reconstruction " +
           "(<b>Analysis&gt;Family history by Dollo parsimony</b>), and " +
           "double-clicking on a family with presence at the ancestral " +
           "node <var>o</var> you are interested in. " +
           "The popup selection menu includes the option " +
           "<var>o</var>&ge; 1 " +
           "if you clicked on a cell with " +
           "value '1' in the column <var>o</var>. " +
           "Extract the filtered rows into a new table " +
           "(<b>Data&gt;Extract selected families...</b>), " +
           "and calculate PGL on that table only.</p>" +
           "<p>The tree display in PGL has an overlaid legend for edge " +
           "length. The legend can be disabled by deselecting the " +
           "<b>Legend</b> checkbox " +
           "in the bottom tool bar.</p>" +
           "<p>Analysis results can be exported into a TAB-delimited text " +
           "file through the popup menu for the corresponding node in " +
           "the Data browser.</p>";
   public static final String CHAPTER_REFERENCES_T = "References";
   public static final String METHOD_REFERENCES 
           =
           "<li>Cs&#369;r&ouml;s, M. (2008). " +
           "\"Ancestral reconstruction by asymmetric Wagner parsimony over " +
           "continuous characters and squared parsimony over distributions.\" " +
           "<em>Springer Lecture Notes in Bioinformatics</em> " +
           "<b>5267</b>:72-86. (Proceedings of RECOMB Satellite Workshop on Comparative Genomics.)</li>" +
           "<li>Cs&#369;r&ouml;s, M. &amp; I. Mikl&oacute;s (2006)." +
           "\"A probabilistic model for gene content evolution with duplication, " +
           "loss, and horizontal transfer.\" " +
           "<em>Springer Lecture Notes in Bioinformatics</em> " +
           "<b>3909</b>:206-220. (Proceedings of RECOMB.)</li>" +
           "<li>Cs&#369;r&ouml;s, M. &amp; I. Mikl&oacute;s (2009a). " +
           "\"Mathematical framework for phylogenetic birth-and-death models.\" " +
           "Report q-bio.PE/0902.0970 at arXiv. " +
           "<a href=\"http://arxiv.org/abs/0902.0970v1\">http://arxiv.org/abs/0902.0970v1</a>.</li>" +
           "<li>Cs&#369;r&ouml;s, M. &amp; I. Mikl&oacute;s (2009b). " +
           "\"Streamlining and large ancestral " +
           "genomes in Archaea inferred with a phylogenetic birth-and-death model.\" " +
           "<em>Molecular Biology and Evolution</em>, <b>26</b>: 2087--2095. " +
           "DOI: 10.1093/molbev/msp123</li>" +
           "";
   public static final String CHAPTER_REFERENCES 
           = "<ul>" +
           METHOD_REFERENCES +
           "<li>Farris, J. S. (1970). " +
           "\"Methods for computing Wagner trees.\" " +
           "<em>Systematic Zoology</em>, <b>19</b>:83Ð92.</li>" +
           "<li>Krylov D. M., Y. I. Wolf, I. B. Rogozin, E. V. Koonin (2003). " +
           "\"Gene loss, protein sequence divergence, gene dispensability, " +
           "expression level, and interactivity are correlated in " +
           "eukaryotic evolution.\" " +
           "<em>Genome Research</em>, " +
           "<b>13</b>:2229Ð2235.</li>" +
           "<li>Yang, Z. (1994). " +
           "\"Maximum likelihood phylogenetic estimation from DNA sequences with " +
           "variable rates over sites: approximate methods.\" " +
           "<em>Journal of Molecular Evolution</em> " +
           "<b>39</b>:306Ð314.</li>" +
           "</ul>";
}