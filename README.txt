	NAME
    streamplot

DESCRIPTION
    Name: Streamer
    Author: rotem.tal
    Description: A tool for printing and clustring streamplots

USAGE
	Run this script from the command line with no arguments to create a configuration file via user friendly interface
	
	
	streamplot.py [-h] [--configuration CFG] [--subject_matrix SUB]
                     [--delim DELIM1] [--meta META] [--meta_delimiter DELIM2]
                     [--meta_sheet SHEET] [--cluster C C C C] [--plot P]
                     [--folder F] [--events [EVENTS [EVENTS ...]]]

	optional arguments:
	  -h, --help            show this help message and exit
	  --configuration CFG   Path to StreamPlot Configuration file (.SCG)
	  --subject_matrix SUB  Path to subject microbiome reads
	  --delim DELIM1        Optional: Subject matrix delimiter (default is ",")
	  --meta META           Path to Meta-data file (as csv, xls, xlsx)
	  --meta_delimiter DELIM2
							Metadata delimiter (default is ",")
	  --meta_sheet SHEET    Metadata sheet name (for Excel files)
	  --cluster C C C C     Cluster subjects and plot clusters, parameters: Number
							of dominant bacteria to consider Matching threshold
							Time threshold Delta
	  --plot P              Plot all inserted subjects
	  --folder F            Folder to save to, default is "Plots"
	  --events [EVENTS [EVENTS ...]]
							Add events from the meta data. Optional arguments (as
							ordered): Event name (default "Onset") Meta data
							subject ID column (default "record_id") Event columns,
							Time column and events identifier (e.g. "Symptomatic",
							"Case" etc.)

CLASSES
    builtins.object
        StreamPlot

    class StreamPlot(builtins.object)
     |  Methods defined here:
     |
     |  __init__(self, subject_name_format='sub_$R_$N_$T', subject_name_sep='_', cmap=None, consts=None, print_other=False, other_n='Other', other_color='k', bacteria_taxo_foramt='$S$B__', bacteria_taxonomic_sep='|', depth='p', depth_dict=None, rec_col='record_id')
     |      Initilize a stream plot object
     |      :param subject_name_format: Format of samples name, $R denotes the ID, $N the sample number, $T sample classification (time stamp, sick visits etc..)
     |      :param subject_name_sep: Separator used in sample names
     |      :param cmap: Color map (as dictionary)
     |      :param consts: Mapping from sample number to time as dictionary (e.g. sample 4 is 6 months, than the entry '4': 6 should be inserted)
     |      :param print_other: Print bacteria not present in the color map
     |      :param other_n: Name for bacteria not present in the color map
     |      :param other_color: Color for bacteria not present in the color map
     |      :param bacteria_taxo_foramt: Bacteria name format, %S denots separator, $B denotes taxa name
     |      :param bacteria_taxonomic_sep: Separator between bacteria taxa ranking
     |      :param depth: Depth at which gradient would be printed
     |      :param depth_dict: Dictionary for possible taxonomic rankings (default is from kingdom to species)
     |      :param rec_col: Meta-Data record ID column name
     |
     |  add_event(self, date, txt, sub_name=None)
     |      Add an event to be printed on a subjects plot
     |      :param date: Location on the stream graph
     |      :param txt: Text for the event (short)
     |      :param sub_name: Subject to add event to
     |      :return:
     |
     |  add_meta(self, meta_csv='data/GMAPMetaClean.xlsx', delim=',', sheet_name='Metadata', header=0, index=0)
     |      Add meta data
     |      :param meta_csv: Path to csv, xls or xlsx
     |      :param delim: data delimiter
     |      :param sheet_name: for xlsx files
     |      :param header: Header row
     |      :param index: Index column
     |      :return:
     |
     |  binary_ordered_cluster(self, k, thr, t_thr=None, clean=False, folder='Plotter')
     |      Clusters subjects, ignoring the abundance level
     |      :param k: Top bacteria
     |      :param thr: Number of corresponding bacteria per time point
     |      :param t_thr: Number of time points for a subject to be clustered
     |      :param clean: Ignore bacteria not in the cmap
     |      :param folder: Folder path to save to
     |      :return:
     |
     |  cluster(self, k, thr, t_thr, delta, folder)
     |      Clusters subjects according to graphical clustering algorithm
     |      :param k: Top bacteria
     |      :param thr: Number of corresponding bacteria per time point
     |      :param t_thr: Number of time points for a subject to be clustered
     |      :param clean: Ignore bacteria not in the cmap
     |      :param delta: allowed abundance level differencefor bacteria
     |      :param folder: Folder path to save to
     |
     |  color_by_abundance(self)
     |      Defines the color mapping of the precision depth based on the overall abundance of the taxa in the data
     |      :return:
     |
     |  events_from_meta(self, event_name='Onset', rec_col='record_id', event_col='time_sx', time_col=None, event_ide=None, event_thr=None, meta_csv='data/GMAPMetaClean.xlsx', delim=',', sheet_name='Metadata', header=0, index=0)
     |      Automatically add events from the meta data
     |      :param event_name: General name for the event (to be printed)
     |      :param rec_col: Meta data ID column name
     |      :param event_col: Event column name
     |      :param time_col: Time column name (to extract event from)
     |      :param event_ide: Match to indicate event
     |      :param event_thr: Threshold (for numerical events) to be considered an event
     |      :param meta_csv: Path to meta data csv
     |      :param delim: meta data delimiter
     |      :param sheet_name: sheet name for xls files
     |      :param header: meta data header row
     |      :param index: meta data index column
     |      :return:
     |
     |  extract_from_name_format(self)
     |      finds and stores the indices (denoted by the subject name separator) of ID ($R), Number ($N) and tag ($T) in the name format
     |
     |  insert_subject(self, sub_csv, header=0, index_col=0, delim=',')
     |      Insert a single subject from csv file
     |      :param sub_csv: Path to csv
     |      :param header: Header row
     |      :param index_col: Index column
     |      :param delim: Delimiter (for csv)
     |      :return:
     |
     |  load_subjects(self, sub_csv, header=0, index_col=0, delim=',')
     |      Insert a csv of several subjects
     |      :param sub_csv: Path to csv
     |      :param header: Header row
     |      :param index_col: Index column
     |      :param delim: csv delimiter
     |      :return:
     |
     |  plot_stream(self, sub_names=None, folder=None, title_format='sub %s streamplot', plot_type='wiggle', suffix='png', x_lab='Time (m)', y_lab='Relative abundance', event_color='#380814')
     |      Plot streamsplots
     |      :param sub_names: List of subjects to plot
     |      :param folder: Folder to save the plots
     |      :param title_format: Plot title format (%s where the subject ID should be inserted)
     |      :param plot_type: Plot baseline style
     |      :param suffix: suffix to save the plots by (file type)
     |      :param x_lab: Label for X axis
     |      :param y_lab: Label for y axis
     |      :param event_color: Color to print events at
     |      :return:
     |
     |  remove_subject(self, sub_name)
     |      Delete a certain subject from streamer memory
     |      :param sub_name: subject name
     |      :return:
     |
     |  save(self, path)
     |      Save current stackplot object
     |      :param path: Path to save to
     |      :return:
     |
     |  set_other(self, print_other, n, c)
     |      Change the Others (bacteria not in the cmap) settings
     |      :param print_other: Boolean, indicate to print others or not
     |      :param n: Others printed name
     |      :param c: Color
     |      :return:
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  DEF_SUB_N = 'DEFAULT_SUB_NM'
     |
     |  UNOCCURRING_STR = 'NOTHINGNEW'

FUNCTIONS
    construct_streamer(streamer)
        Construct the streamer object
        :param streamer: object to be constructed or modified
        :return:

    get_name(streamer)
        Get subject name format from the user
        :param streamer: Current streamer object
        :return:

    main()
        Run this program without arguments to configure a new streamer object
        :return:

    read_dict(some_dict)
        Reads a colormap dictionary from json, csv or a string representing a
        dictionary
        :param some_dict: Path to .json or .csv file, or string
        :return:

    run(streamer)
        Run the configuration editor
        :param streamer: Streamer object to modify
        :return:

DATA
    rcParams = RcParams({'_internal.classic_mode': False,
         ...nor.widt...
    rcParamsDefault = RcParams({'_internal.classic_mode': False,
         ...n...
    rcParamsOrig = RcParams({'_internal.classic_mode': False,
         ...nor....


