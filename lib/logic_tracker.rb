# A class that allows for extensible logic tracking of decisions.
#
# The intended usage of this class is to create an object that 
#   can incrementally generate unique identifiers as keys
#   with values being a text string that the programmer can define.
#   The final output is to concatenate all of the identifiers into
#   into a delimited string (user defined at instantiation) to be
#   placed into a VCF INFO tag, and a key-value definition text file
#   externally written for decoding and review of the delimited string.
#
# Goal of LogicTracker
#   The goal here is to allow the developer to place n number of unique trackers with
#   descriptions without needed to keep track of them all themselves. The unique names
#   may be as verbose as the developer requires without complicating compact output and
#   recording of trackers added to any instance of LogicTracker.
#
# Terminology:
#   track/tracker: an instance of a Trail class that records that a unique Trail was passed/added.
#   trail: a series of tracks that makeup the history of markers/tracks
#   Trail: an internal use object only that retains the name, id, and description of a track.
#     Each Trail instance should be unique and the ids should be consecutive unless an error occurred.
#
# Using logic_tracker and structure
#   The LogicTracker class retains all unique Trail instances at a class level. This ensures
#   that each Trail is unique (by both name and id, separately) so that any Trail name or
#   id may implicates only 1 decision, ever. Also at a class level, is an internally incrementing i
#   variable used to identify each Trail instance without the user needing to be concerned.
#   This class level internal i variable increments when a new Trail instance (aka. track) is 
#   created via the add_track function. Also, a broken consecutive series of id
#   values would indicate that a name had somehow been overwritten and that name-id-description
#   is not unique. The 
#
#   At a LogicTracker instance level, a 'trail' of tracks (aka. Trail instances) added is retained.
#   This trail represents all of the Trail markers passed in recent LogicTracker instance history.
#   At this level, there is an add_track function that will append a track to the instance's record
#   and, if that track's name has never been seen before in the LogicTracker's class, it will be
#   recorded as part of all possible Trails. Finally, there are various output functions but most
#   important is get_trail function. This function allows outputs the instance's Trail history in
#   either a detailed or compacted format. The compacted format is a string of Trail ids, in order
#   of being add_tracked, delimited by a user defined delimiter when instantiating that particular 
#   instance of LogicTracker. The detailed option outputs similar information but expands on it by
#   providing the Trail markers/tracks, in order of being add_tracked, with the Trail's name, id
#   and description.
#
# Analogy To Help Explain:
#   A mountain pass (the program) has many paths through it. Some are traversed and marked
#   for every and all hikers (LogicTracker class level). Each individual hiker 
#   (LogicTracker instance) may pass some trail markers (Trails) while placing new ones for new paths. 
#   A record is of all markers any one hiker may see on their hike is kept for use by all hikers.
#   Each individual hiker records all trail markers they have seen, and add new markers when they 
#   pathfind toward never before traversed trails.
#
# Example: Kafeen classification tracking
#   An instance of the Command class retains an instance of LogicTracker. Each pathogenicty
#   classification determining decisions has an add_track inserted into it so that, the 1st
#   time that decision is made it is added to the class level list of all Trails, given a unique
#   id and then added to the LogicTracker's instance list of Trails seen. Any previously
#   seen/known Trails that are add_tracked, are only added to the LogicTracker's instance list.
#   Both add_predictions and finalize_predictions have this tracker integrated into only the
#   pertinent parts of its logic. In the add_predictions function, each variant assessed
#   eventually has its compact list of Trails inserted into an INFO tag (CLASSIFY_TRAIL_PART)
#   via use of the get_trail function. The LogicTracker instance's trail is cleared at the end 
#   of each variant's assessment to provide a 'clean' trail for the next variant. The 
#   finalize_pathogenicity function, does something similar in that each variant is assessed 
#   and output to an INFO tag (CLASSIFY_TRAIL). However, finalize_pathogenicity reads in the 
#   CLASSIFY_TRAIL_PART for that variant and appends its add_track Trails to it for its final 
#   output. At the very end of the finalize_pathogenicity function, the entire list of Trails 
#   that have been seen during that run of Kafeen is output to FILE_PREFIX.classify_track_list.tsv. 
#   By using the FILE_PREFIX this file keeps with naming conventions and is written to the directory 
#   that all other Kafeen files using the FILE_PREFIX are.
#
# @author Rob Marini

require 'logger'
require 'fileutils'
require_relative 'core_extensions'

# Monkey-patch the String class
String.include CoreExtensions::String::Vcf
String.include CoreExtensions::String::Colorize

class LogicTracker
  
  #class variables
  @@i = 1;  # auto_incrementing index, user begins at 1, 0 is reserved for non-decisions/null decision/no decision made.
  @@trackers = Hash.new; # the key value pair tracker, key is a symbol name, value is an array of [id,description]
  
  #class functions
  def get_all_trackers()
    return(@@trackers)
  end
  
  def list_trackers()
    tracker_list = [["NAME\tID\tDESCRIPTION\n"]]
    for tracker in @@trackers do
      tracker_list.push(tracker[1].list_members("\t") + "\n")
    end
    return(tracker_list)
  end
  
  def get_i_incr()
    i = @@i
    @@i += 1
    return i
  end
  
  def initialize(delim, tracking:true, log_level:'info', log_out:STDOUT)
    
    # Set Logger
    @@log = Logger.new(log_out)
    if log_level.upcase == 'UNKNOWN'
      @@log.level = Logger::UNKNOWN
    elsif log_level.upcase == 'FATAL'
      @@log.level = Logger::FATAL
    elsif log_level.upcase == 'ERROR'
      @@log.level = Logger::ERROR
    elsif log_level.upcase == 'WARN'
      @@log.level = Logger::WARN
    elsif log_level.upcase == 'INFO'
      @@log.level = Logger::INFO
    elsif log_level.upcase == 'DEBUG'
      @@log.level = Logger::DEBUG
      @@log.debug("Debugging messages enabled")
    else
      @@log.level = Logger::DEBUG
      @@log.error("#{log_level} is not a valid log level")
      @@log.warn("Debugging messages enabled by default")
    end
    
    #object variables
    @delim = delim  # delimiter for concatenation of indexes
    @trail = Array.new; #empty array for retaining the decisions made throughout the run. Only retains the id's of trackers
    @@trackers["void".to_sym] = Trail.new("void", 0, "No decision was made here. Consider this like a null decision.")
    
  end
  
  #instance functions
  def add_track(logicName, note='NO DESCRIPTION PROVIDED')
    if(@@trackers.key?(logicName.to_sym))
      #key exists, just add the track
      @trail.push(@@trackers[logicName.to_sym].get_id)
    else
      #key doesn't exist, add it and the track
      @@trackers[logicName.to_sym] = Trail.new(logicName, self.get_i_incr(), note)
      @trail.push(@@trackers[logicName.to_sym].get_id)
    end
  end
  
  def clear_trail()
    @trail = Array.new;
  end
  
  def get_trail(detail=false)
    if(detail)
      trail_str = [["Logical decisions taken in order:"],["Compact form: " + @trail.join(@delim)],["NAME","ID","DESCRIPTION"].join("\t")]
      for id in @trail do
        trail_str.push(@@trackers[@@trackers.keys[id]].list_members("\t"))
      end
      return(trail_str.join("\n"))
    else
      return(@trail.join(@delim))
    end
  end
  
end

# internal-only simple class for retaining name, id, and description information
class Trail
  def initialize(logicName, id, desc)
    @name = logicName
    @id = id
    @desc = desc
  end
  
  def get_id()
    return @id
  end
  
  def get_name()
    return @name
  end
  
  def get_desc()
    return @desc
  end
  
  def list_members(delim)
    return([@name,@id,@desc].join(delim))
  end
  
end