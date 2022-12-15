import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


import re
import copy
import math
import numpy as np
import pandas as pd

import baltic3 as bt

# date conversion imports
import datetime
from dateutil.relativedelta import relativedelta as rd
from datetime import timedelta
import time

"""A bunch of functions which I wrote to support baltic3.py.
"""

def get_leaf(tree, tipname):
    """Gets a leaf from `tree` by its tipname.

    PARAMS
    ------
    tree: baltic tree object.
    tipname: str; name of tip.

    RETURNS
    -------
    lf: baltic leaf object with name `tipname`
    """
    lf = None
    for k in tree.leaves:
        if k.name == tipname:
            lf = k
    return lf


def quick_draw_tree(tree, 
                    colour_by = "",
                    values_of_interest = [],
                    fig_h = 12,
                    fig_w = 9,
                    branch_width=0.5, 
                    branch_colour="gray",
                    tip_shape_size=8, 
                    x_offset=0, 
                    save_fn = "",
                    verbose=True,
                    show_borders=False):
    """
    Draws a tree, and colours a set of tips based on an (optional) input dataframe, dm, with a column of interest, `colname`. 
    
    Params
    ------
    tree: input baltic tree.
    colour_by: str; trait key name to colour the tips by. Can be left blank.
    values_of_interest: list of str; optional. a subset of all possible values of the trait of interest, 
    e.g. only "human" and "avian" out of all possible hosts in the dataset, or only "USA" and "China" 
    out of all possible countries, etc. 
    If `colour_by` is an empty string, this parameter does nothing, even if specified with a non-empty list.
    If left as an empty list, will grab all possible values in the `colour_by` trait column value. 
    WARNING: The latter option may result in too many possible colour values to map. 
    fig_h: float; figure height.
    fig_w: float; figure width.
    branch_width: float; branch width (line weight).
    branch_colour: str; branch colour. Default black. 
    tip_shape_size: float; tip shape size.
    x_offset: float; x_offset. The root is placed at x=0 by default. 
    save_fn: str; output filename. Does not save if left blank.
    verbose: Boolean; verbosity param
    show_borders: Boolean; if true, show axes and border elements. 

    Returns
    -------
    Tree plot on the active notebook. 
    """

    # Init c_dict
    c_dict = {}
    if colour_by != "":
	    if len(values_of_interest) > 0:
	        c_dict = brew_colour_dictionary(values_of_interest)
	    else:
	        all_trait_vals = [lf.traits[colour_by] for lf in tree.leaves]
	        c_dict = brew_colour_dictionary(all_trait_vals)
        
    # ==================== Plot! ====================
    fig,ax = plt.subplots(figsize=(fig_w, fig_h),facecolor='w')

    for k in tree.Objects:
        x=k.height
        y=k.y

        xp = k.parent.height
        if x is None:
            x = x_offset
        if xp==None:
            xp = x + x_offset

        if isinstance(k,bt.leaf) or k.branchType=='leaf':
            if colour_by != "":
                # colour only those tips identified in c_dict
                if k.traits[colour_by] in list(c_dict.keys()):
                    ax.scatter(x, y, facecolor=c_dict[k.traits[colour_by]], 
                               edgecolor='none',
                               s=tip_shape_size, 
                               zorder=10)
                    # following Gytis' art style: plot black circles underneath
                    ax.scatter(x,y,s=tip_shape_size+0.8*tip_shape_size,
                               facecolor='k',
                               edgecolor='none',
                               zorder=9)
            else:
                pass

        elif isinstance(k,bt.node) or k.branchType=='node':
            ax.plot([x,x],
                    [k.children[-1].y,k.children[0].y],
                    lw=branch_width,
                    color=branch_colour,ls='-',zorder=9)

        # Draw horizontal lines
        ax.plot([xp,x],[y,y],lw=branch_width,color=branch_colour,ls='-',zorder=9)

    # ==================== Figure Legend ====================
    if colour_by != "":
        labels = [] # for ax.legend()
        for c in list(c_dict.keys()):
            labels.append(mpatches.Patch(color=c_dict[c], label=c))

    ax.legend(handles=labels)

    # scale bar
    scalebar_y = -tree.ySpan*0.05
    ax.plot([0, 0.01], [scalebar_y, scalebar_y], c="k", lw=branch_width*2)
    ax.text(0.005, -tree.ySpan*0.045, "0.01", 
            verticalalignment="bottom", 
            horizontalalignment="center")
    
    # ==================== remove tick marks and borders ==================
    if "show_borders"==False:
        ax.set_yticks([])
        ax.set_xticks([])
        plt.axis('off')
    #plt.tight_layout()    

    # Save the figure to a pdf. 
    if save_fn != "":
        plt.savefig(save_fn, bbox_inches="tight")

    plt.show()


def austechia_read_tree(tree_path, date_bool=False, date_pos=-1, date_delim="_", make_tree_verbose=False):
    """Lifted from the austechia.ipynb (thus the name). 
    This works for BEAST and RAXML trees, or raw newick strings, but not really for treetime or LSD dated trees. 
    This is more because of the output: treetime-dated trees don't have the absoluteTime directly written in the newick strings. 
    In those cases, manually extrapolate the node times from known tip dates. The branch lengths, at least, will be proportionate
    as calendar times. 
    
    A number of lines which parsed dates in each tipname commented out,
    So now absoluteHeight attributes should be empty.
    Tested and working on: sample treesub output, sample mcc tree.

    PARAMS
    ------
    tree_path: str; path to input free file.
    date_bool: Bool; checks if dates are in the tipnames. If False, `date_pos` and `date_delimiter` input arguments will do nothing.
    date_pos: int; index of the dates in tipnames, if available (specify with `date_bool`). Decimal or calendar date formats are both supported. 
    Calendar dates must be in yyyy-mm-dd, yyyy-mm or yyyy format. 
    date_delimiter: delimiter vale to read date off the tipname. 
    make_tree_verbose: Bool; verbosity parameter.

    RETURNS
    -------
    ll: baltic tree object. 
    """
    tipFlag=False
    tips={}

    for line in open(tree_path,'r'): ## iterate through FigTree lines
        l=line.strip('\n') ## strip newline characters from each line

        cerberus=re.search('dimensions ntax=([0-9]+);',l.lower()) ## check how many tips there are supposed to be
        if cerberus is not None:
            tipNum=int(cerberus.group(1))

        #####################
        cerberus=re.search('tree TREE([0-9]+) = \[&R\]',l) ## search for beginning of tree string in BEAST format
        if cerberus is not None:
            treeString_start=l.index('(') ## tree string starts where the first '(' is in the line
            ll=bt.tree() ## new instance of tree
            bt.make_tree(l[treeString_start:],ll, verbose=make_tree_verbose) ## send tree string to make_tree function, provide an empty tree object
        #####################

        if tipFlag==True:
            cerberus=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l) ## look for tip name map, where each tip is given an integer to represent it in tree
            if cerberus is not None:
                tips[cerberus.group(1)]=cerberus.group(2).strip("'") ## if you give tips an integer (in the form of a string), it will return the full name of the tip
            elif ';' not in l: ## something's wrong - nothing that matches the tip regex is being captured where it should be in the file
                print('tip not captured by regex:',l.replace('\t',''))

        if 'translate' in l.lower(): ## start looking for tips
            tipFlag=True
        if ';' in l: ## stop looking for tips
            tipFlag=False

    #print("Number of objects found in tree string: %d"%(len(ll.Objects)))

    ## rename tips, find the highest tip (in absolute time) in the tree
    if len(tips)==0:
        for k in ll.Objects:
            if isinstance(k,bt.leaf):
                k.name=k.numName

        # read the tip date. Accepts decimal or calendar dates;
        # converts to a decidate if required.
        if date_bool:
            highestTip=max([decimalDate(x.name.strip("'").split(date_delim)[date_pos],variable=True) for x in ll.Objects if isinstance(x,bt.leaf)])
    else: ## there's a tip name map at the beginning, so translate the names
        ll.renameTips(tips) ## give each tip a name
        if date_bool:
            highestTip=max([decimalDate(x.strip("'").split(date_delim)[date_pos],variable=True) for x in tips.values()])

    ll.treeStats()
    ll.sortBranches(descending=False)
    if date_bool:
        ll.setAbsoluteTime(highestTip)
        print('Highest tip date: %.4f'%(highestTip))


    return ll


def loadNexus(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',treestring_regex='tree [A-Za-z\_]+([0-9]+)',variableDate=True,absoluteTime=True,verbose=False):
    """Gytis' original tree-reading function.
    """
    tipFlag=False
    tips={}
    tipNum=0
    ll=None
    if isinstance(tree_path,str):
        handle=open(tree_path,'r')
    else:
        handle=tree_path

    for line in handle:
        l=line.strip('\n')

        cerberus=re.search('dimensions ntax=([0-9]+);',l.lower())
        if cerberus is not None:
            tipNum=int(cerberus.group(1))
            if verbose==True:
                print('File should contain %d taxa'%(tipNum))

        cerberus=re.search(treestring_regex,l)
        if cerberus is not None:
            treeString_start=l.index('(')
            ll=bt.tree() ## new instance of tree
            bt.make_tree(l[treeString_start:],ll) ## send tree string to make_tree function
            if verbose==True:
                print('Identified tree string')

        if tipFlag==True:
            cerberus=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l)
            if cerberus is not None:
                tips[cerberus.group(1)]=cerberus.group(2).strip('"').strip("'")
                if verbose==True:
                    print('Identified tip translation %s: %s'%(cerberus.group(1),tips[cerberus.group(1)]))
            elif ';' not in l:
                print('tip not captured by regex:',l.replace('\t',''))

        if 'translate' in l.lower():
            tipFlag=True
        if ';' in l:
            tipFlag=False

    assert ll,'Regular expression failed to find tree string'
    ll.traverse_tree() ## traverse tree
    ll.sortBranches() ## traverses tree, sorts branches, draws tree
    if verbose:
        print("len(tips)=%s" % len(tips))
    if len(tips)>0:
        ll.renameTips(tips) ## renames tips from numbers to actual names
        ll.tipMap=tips
    if absoluteTime==True:
        tipDates=[]
        for k in ll.leaves:
            if len(tips)>0:
                n=k.name
            else:
                n=k.numName

            cerberus=re.search(tip_regex,n)
            if cerberus is not None:
                tipDates.append(decimalDate(cerberus.group(1),fmt=date_fmt,variable=variableDate))

        highestTip=max(tipDates)
        ll.setAbsoluteTime(highestTip)

    return ll


def treesub_to_bt(fn_in, fn_out, verbose=True):
    """
    IMPT NOTE: dm output not working. Parse substitutions.tsv output directly instead
    for the moment...
    Reads raw treesub output, which is assumed to be in nexus format, does a
    bunch of necessary wrangling, then writes out a format that can be read by
    austechia_read_tree(). 

    The internal nodes of raw treesub output tree have too many parameters that are 
    poorly delimited, which makes it difficult for baltic's regex to read. We only need 
    two parameters anyway: the unique node identifier (NUMBER) and the NS substitutions 
    (nonsynsub). The new tree written out to file will only retain NUMBER; the NS subs
    can be pulled from the daataframe generated by this function. 

    Params
    ------
    fn_in: string; treesub output, like "substitutions.tree"
    fn_out: filename to write the modified nexus format out to.

    Returns
    -------
    dm: dataframe of node_num (originally "NUMBER") and the associated non-syn
    subs.

    More Notes
    ----------
     * I may have overdone the "raw treesub output" bit; this function can't
     read trees rerooted in figtree (but I haven't tried very hard either).
     Treesub must be rerun with a new root.
    """
    # Assumes that raw treesub output is in nexus format
    # read treestring
    with open(fn_in) as f:
        contents = f.readlines()
    contents = [x.strip() for x in contents]

    # Get the row which contains just the newick string
    for i in range(len(contents)):
        if contents[i] == "begin trees;":
            tree_string = contents[i+1]

    # Remove all spaces, and a bit of nexus formatting in front
    tree_string = tree_string.replace(" ", "")
    tree_string = tree_string.replace("treetree_1=[&R]", "")

    # =============== Extract stuff ===============
    # Get tip names
    tip_names_ls = []
    tip_name_start = 0
    for i in range(len(contents)):
        # switch on
        if contents[i] == "taxlabels":
            tip_name_start = 1
        # switch off
        if (contents[i] == ";") and (contents[i+1] == "end;") and (tip_name_start == 1):
            tip_name_start = 0

        if (tip_name_start == 1) and (contents[i] != "taxlabels"):
            pattern = re.compile(r"([\s\S]*?)\[&")
            tipname = re.findall(pattern, contents[i])[0]
            tip_names_ls.append(tipname)

    # Retrieve all treesub node_strings, "[&REALNAME=...]", from tree_string
    pattern = re.compile(r'\[&([\s\S]*?)\]:')
    node_strings_ls = re.findall(pattern, tree_string)

    # Replace the full treesub nodestring with just "NUMBER"
    # Keep a dictionary
    for node_string in node_strings_ls:
        # retrieve "NUMBER"
        pattern = re.compile(r'NUMBER="[\d]*?"')
        node_num_ls = re.findall(pattern, node_string)
        node_num = ""
        if len(node_num_ls) > 0:
            node_num = node_num_ls[0]
        else:
            print("WARNING: NUMBER not found for node_string %s" % node_string)
        tree_string = tree_string.replace(node_string, node_num)

    # Extract all the individual bits of node_strings
    contents = []
    for node_string in node_strings_ls:
        #old NUMBER regex: r'NUMBER="([\d]*?)",'
        pattern = re.compile(r'NUMBER="[\d]*?"')
        node_num = re.findall(pattern, node_string)[0]

        pattern = re.compile(r'NONSYNSUBS="([\s\S]*?)",')
        nssubs_ls = re.findall(pattern, node_string)
        nssubs = "*"
        if len(nssubs_ls) > 0:
            nssubs = nssubs_ls[0]
            for char in ["[", "]"]:
                nssubs = nssubs.replace(char, "")
        contents.append([node_num, nssubs])

    dm = pd.DataFrame(data=contents, columns=["node_num", "nonsynsubs"])
    # some formatting
    dm["node_num"] = dm.apply(lambda row: str(row["node_num"]).replace('NUMBER="', ""), axis=1)
    dm["node_num"] = dm.apply(lambda row: str(row["node_num"]).replace('"', ""), axis=1)

    # =============== Prep contents to write out a full nexus file ===============
    contents = ["#NEXUS"]
    for ln in ["begin taxa;", "dimensions ntax="+str(len(tip_names_ls))+";", "taxlabels"]:
        contents.append(ln)

    # add taxlabels
    for ln in tip_names_ls:
        contents.append(ln)

    contents.append(";")
    contents.append("end;")

    # add tree block
    tree_block_idx = len(contents)
    contents.append("begin trees;")
    contents.append("tree TREE1 = [&R] "+tree_string)
    contents.append("end;")

    with open(fn_out, "w") as f:
        for line in contents:
            f.write("%s\n" % line)

    # preview contents
    if verbose:
        preview_counter = 10
        for i in range(preview_counter):
            print(contents[i])
        print("...")
        print(contents[tree_block_idx])
        print(contents[tree_block_idx+1][:50]+"...")
        print(contents[tree_block_idx+2])
        print("")
        print("Written out to file %s" % fn_out)

    return dm


def assign_leaf_trait(tree, dm, query_colname, target_colname, trait_name=""):
    """
    Assigns the values of target_colname to each leaf trait dictionary in the input tree.

    Params
    ------
    tree: input Baltic tree
    dm: pandas dataframe
    query_colname: str; the column name from which to look up each leaf name.
    target_colname: str; the column name of the trait of interest.
    trait_name: str; the name of the dictionary key. If left blank (default), will inherit the value of `target_colname`.

    Returns
    -------
    tree: tree with leaf traits assigned in-place.
    """
    if trait_name == "":
        trait_name = target_colname

    trait_val = "undef"
    for lf in tree.leaves:
        d_t = dm.loc[dm[query_colname]==lf.name]
        if len(d_t) != 1:
            print("WARNING: %s records found in dataframe for leaf name %s! Will assign 'undef'." % (len(d_t), lf.name))
        elif len(d_t) == 1:
            trait_val = list(d_t[target_colname])[0]
        lf.traits[trait_name] = trait_val

    return tree


def assign_inode_traits(tree, trait_name):
    """
    Iterate over all internal nodes, assigning traits to each inode based on their leaves.
    If all leaves of that node have the same trait, then that inode will have that trait.
    Otherwise, that node will be assigned 'undef'. Helpful in colouring the branches of monophyletic clades:
    branch colours are inherited from the trait colour assignment of each parent node. 
    This code is actually pretty inefficient, because the same node.leaves operation is called for _all_
    inodes. 
    
    Params
    ------
    tree: Baltic tree with some trait of interest assigned to all leaves.
    trait_name: str; trait name to look up in each leaf.traits dictionary. Will be assigned as a key in each node.traits
    dictionary.

    Returns
    -------
    tre2: Baltic tree with nodes assigned with traits. 
    """

    inodes_ls = tree.nodes

    inode_trait = "undef" # init
    for nd in inodes_ls:
        leaves_ls_temp = nd.leaves # Leaves of the current node
        lf_traits_ls = []
        for lf in nd.leaves:
            lf_traits_ls.append(lf.traits[trait_name])

        if len(np.unique(lf_traits_ls)) == 1:
            node_trait = np.unique(lf_traits_ls)[0]
        
        # write nd.traits in-place, which I hate.
        nd.traits[trait_name] = node_trait

    return tree


def brew_colour_dictionary(my_list, scheme="qualitative", style="paired"):
    """Create a colour dictionary based on colorbrew presets.
    
    Params
    ------
    my_list: a list of values that will be assigned one colour each.
    scheme: str; sequential, diverging or qualitative. Dummy input that does nothing (yet).
    style: str; only 'paired' and 'set1' from 'qualitative' scheme available. 
    """
    cbrew_full_dict = {"paired":["#a6cee3", "#1f78b4", "#b2df8a",
                                 "#33a02c",
                                 "#fb9a99",
                                 "#e31a1c",
                                 "#fdbf6f",
                                 "#ff7f00",
                                 "#cab2d6",
                                 "#6a3d9a",
                                 "#ffff99",
                                 "#b15928"],
                       "set1":["#e41a1c", "#377eb8", "#4daf4a", 
                               "#984ea3", #purple
                               "#ff7f00", #orange
                               "#ffff33", #yellow
                               "#a65628", #brown
                               "#f781bf", #pink
                               "#999999"] #grey
                      }
    
    full_colour_list = cbrew_full_dict[style]
    cdict = {}
    try:
        for idx, val in enumerate(my_list):
            cdict[val] = full_colour_list[idx]
            
    except IndexError:
        print("Too many elements in input list to map the colour style to!")
        cdict = {}
        
    return cdict

def preorder_traverse(node, preorder_ls, verbose=True):
    """Recursive algorith to traverse a tree, given a starting (current) node.
    Typical usage: input my_tree.cur_node, where my_tree is a baltic tree object.
    This tells preorder_traverse() to start traversal from the root.
    Prints out the order in which the nodes are visited.

    Usage:
    >>> preorder_ls = [] # init list to fill
    >>> preorder_traverse(my_input_node, preorder_ls)
    """
    preorder_ls.append(node)
    if node.branchType == "node":
        if verbose:
            print(node)
        preorder_traverse(node.children[0], preorder_ls)
        preorder_traverse(node.children[1], preorder_ls)
    elif node.branchType == "leaf":
        if verbose:
            print(node, node.name)


def postorder_traverse(node):
    opers = {'+':operator.add, '-':operator.sub, '*':operator.mul, '/':operator.truediv}
    res1 = None
    res2 = None
    if tree:
        res1 = postordereval(node.children[0])
        res2 = postordereval(node.children[1])
        if res1 and res2:
            return opers[node.parent](res1,res2)
        else:
            return node.parent


def decimalDate(date,fmt="%Y-%m-%d",variable=False,dateSplitter='-'):
    """ Converts calendar dates in specified format to decimal date.

    If the input date is already a float, does nothing.
    THIS IS NOT AN ADVISABLE HACK BECAUSE IT HIDES ERRORS.
    """

    if isfloat(date):
        #print("date already a decimal year, pass.")
        result = float(date)
    else:
        if variable==True: ## if date is variable - extract what is available
            dateL=len(date.split(dateSplitter))
            if dateL==2:
                fmt=dateSplitter.join(fmt.split(dateSplitter)[:-1])
            elif dateL==1:
                fmt=dateSplitter.join(fmt.split(dateSplitter)[:-2])

        adatetime=datetime.datetime.strptime(date,fmt) ## convert to datetime object
        year = adatetime.year ## get year
        boy = datetime.datetime(year, 1, 1) ## get beginning of the year
        eoy = datetime.datetime(year + 1, 1, 1) ## get beginning of next year

        ## return fractional year
        result = year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds()))

    return result


def isfloat(value):
    """Checks if a string can be converted into a float.
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def datetime_to_decimal_year(date):
    """Some SO solution I haven't really looked into yet.

    PARAMS
    ------
    date: datetime object

    returns
    -------
    decimal year: float. the decimal year version of the input.

    e.g.
    >>> T = datetime.date(2016, 11, 29)
    >>> most_recent = to_decimal_year(T)
    >>> most_recent
    2016.9098360655737
    """
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime.datetime(year=year, month=1, day=1)
    startOfNextYear = datetime.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def to_decimal_year(date_str, date_delimiter="-"):
    """Wraps datetime_to_decimal_year, with the additional precursor step
    of converting the input date_str into a datetime.date object.

    Params
    -------
    date_str: str, date in the format yyyy-mm-dd, yyyy-mm, or yyyy
    date_delimiter: str, separator for the input date_str

    Returns
    -------
    decimal year: float. the decimal year version of the input.
    """
    date_ln = date_str.split(date_delimiter)
    T = "*"
    if len(date_ln) == 3:
        T = datetime.date(int(date_ln[0]), int(date_ln[1]), int(date_ln[2]))
    elif len(date_ln) == 2:
        T = datetime.date(int(date_ln[0]), int(date_ln[1]), 1)
    elif len(date_ln) == 1:
        T = datetime.date(int(date_ln[0]), 1, 1)

    return datetime_to_decimal_year(T)


#def unique(o, idfun=repr):
#    """Reduce a list down to its unique elements."""
#    seen = {}
#    return [seen.setdefault(idfun(e),e) for e in o if idfun(e) not in seen]
