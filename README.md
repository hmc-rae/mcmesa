# mcmesa

An expansion wrapper of the python mesa_reader library to allow easier parsing of MESA stellar model results. Soon to be expanded into its own reader library entirely.
Made by Harris C. McRae. Use at your own risk - early work in progress...

# Requirements

- numpy 2.4.0
- mesa-reader 0.3.5

At a later version of McMesa, the requirement for mesa-reader will be phased out in favor of custom file reading.

# Features

McMesa defines two key objects - **MesaLogs** and **ProfileWrapper**.

**MesaLogs** defines an entire collection of MESA logs, consisting of a history log (with details for all models in the MESA simulation) and multiple profiles, each of which contain stellar zone data for a particular model. A **MesaLogs** object provides direct access to the history log, and also links to loaded profiles for specific models.

**ProfileWrapper** defines in a similar manner the parameters of a specific profile. Multiple **ProfileWrapper** objects are loaded when a **MesaLogs** object is created.

Both **MesaLogs** and **ProfileWrapper** objects provide the same access method to the parameters represented with Python list-access notation. See *Usage* for examples on how to access columns.

A list of features is as follows:
 - Automatic loading & association of profiles with models with the **MesaLogs** object.
 - Optional as-needed loading of profiles within the **MesaLogs** object (only loads from disk when a profile is referenced).
 - Iteration through profiles and profile array slicing on model numbers.
 - Reading of logarithmic columns as both logarithm and linear data at runtime.
 - Methods to facilitate expansion of data in profiles (i.e. adding new columns).
   - Methods to allow constructive expansion of data based on order of some other column (i.e. compute difference in ascending radius).

# Usage

To create a **MesaLogs** object, simply provide the path to the folder containing the following:
 - history.data
 - profiles.index
 - profile(N).data

With an example file system as below:
```
~/MesaModels
    ~/MesaModels/LOGS_to_end_core_h_burn
        ~/MesaModels/LOGS_to_end_core_h_burn/history.data
        ~/MesaModels/LOGS_to_end_core_h_burn/profiles.index
        ~/MesaModels/LOGS_to_end_core_h_burn/profile1.data
        ~/MesaModels/LOGS_to_end_core_h_burn/profile2.data
        ~/MesaModels/LOGS_to_end_core_h_burn/profile3.data
```
loading a **MesaLogs** of this object would be done with the following:
```py
import mcmesa as mcm
CoreHBurn = mcm.MesaLogs('~/MesaModels/LOGS_to_end_core_h_burn')
```

The object `CoreHBurn` now represents on its own the `history.data` file. Data from the `history.data` file can be accessed as shown below:
```py
CoreHBurn.ColumnNames           # Returns a tuple of all column names
CoreHBurn.Shape                 # Returns a tuple for the shape of the history data
CoreHBurn['model_number']       # Returns the model_number column
CoreHBurn['star_age', 50]       # Returns the value of column star_age at model_number=50
CoreHBurn['all', 50]            # Returns a dict object representing all data values at row model_number=50
CoreHBurn[None, 50]             # Same as above
```

The **MesaLogs** object uses the `'model_number'` column as its root index - rows are accessed based on values in the `'model_number'` column unless the attribute `AutoCorrectModelNum` is `False` (default : `True`).

The profiles for this **MesaLogs** object are stored under the attribute `Profiles`, a **MesaLogs.ProfileSet** object. Accessing can be done as such:
```py
CoreHBurn.ModelProfiles             # Returns a tuple of all models that have a profile
CoreHBurn.Profiles.ModelNumbers     # Same function as above - tuple of model numbers for each profile.
CoreHBurn.Profiles[100]             # Returns the profile for model number 100, if it exists.
```

The **MesaLogs.ProfileSet** object under `Profiles` provides access to each profile based on the model number associated. If automatic loading of profiles from disk is disabled (with `autoLoad=False` in the **MesaLogs** constructor), then the profiles will only be loaded once a profile is accessed through the methods described here.

As earlier mentioned, slicing & iterator support is provided for the **MesaLogs.ProfileSet** object:
```py
# Create a new ProfileSet object that has only profiles with model numbers between 100 and 299.
SubProfiles = CoreHBurn.Profiles[100:299]

# Create a new ProfileSet object that only has profiles with model numbers between 100 and 299, selecting every 2nd one.
SubProfiles = CoreHBurn.Profiles[100:299:2]

# Iterates through every profile in the full set of profiles
for profile in CoreHBurn.Profiles:      
    pass

# Iterate through every profile with model numbers between 100 and 299
for profile in CoreHBurn.Profiles[100:299]:      
    pass
```

Accessing elements of the **MesaLogs.ProfileSet** object through either iteration or direct array accessing, as shown above, returns the **ProfileWrapper** object. The **ProfileWrapper** object has core accessing identical to the **MesaLogs** object.

```py

Profile = CoreHBurn.Profiles[100]
Profile.ColumnNames                 # Returns a tuple of all column names
Profile.Shape                       # Returns a tuple for the shape of the history data
Profile['zone']                     # Returns the zone column
Profile['logRho', 50]               # Returns the value of column logRho at zone=50
Profile['Rho']                      # Returns the virtual column Rho (10^logRho)
Profile['all', 50]                  # Returns a dict object representing all data values at row zone=50
Profile[None, 50]                   # Same as above
```

Similar to how the **MesaLogs** object accesses rows by default, the **ProfileWrapper** object accesses rows in this manner by zone number. Zones are indexed from 1. This can also be disabled by setting the **ProfileWrapper** attribute `AutoZone` to be `False` (default : `True`).

Contents of the **ProfileWrapper** data can be assigned to using the same array access notation. Both columns & individual cell values can be set. This functionality is not recommended.

```py
# Fill in the logRho column with all zeroes.
Profile['logRho'] = np.zeros((1, Profile.Shape[1])) 

# Set the value at logRho at zone 1 to be 100.
Profile['logRho', 1] = 0
```

Structured methods to compute new columns exist in the **ProfileWrapper** object, with methods `AddColumn` and `ConstructNewData`.

```py
def double_logRho(zone:dict)->float:
    return zone['logRho'] * 2

ones = np.ones((1, Profile.Shape[1]))
Profile.AddColumn('double_logRho', ones)    # Create a new column in this ProfileWrapper called 'double_logRho' that has the contents of the array 'ones'.

# Create a new column in this ProfileWrapper called 'double_logRho', overwriting the previous one, that has contents defined by the function double_logRho at each zone.
Profile.AddColumn('double_logRho', double_logRho, allowOverwrite=True)
```
If no argument is provided for the column contents, it defaults to zero.

```py
def base_compute(zone):
    return zone['logRho']

def step_compute(current_zone, previous_zone):
    return current_zone['logRho'] + previous_zone['sum_logRho']

# Constructs a new column titled 'sum_logRho' in this ProfileWrapper. 
# This column is the sum of logRho throughout the profile.
# The new column is constructed in ascending order of logR (log Radius) in ascending order.
# The lowest radius zone will be computed initially according to base_function, and every column after that will be computed by step_compute, which takes two zones for the current and previous zone datas.
Profile.ConstructNewData('sum_logRho', 'logR', base_function=base_compute, step_function=step_compute, direction='asc')
```

Finally, **ProfileWrapper** objects can access its associated parent **MesaLogs** object with the attribute `History`, and directly access its row in the **MesaLogs** object with the method `HistoryState(column)`.

```py
Profile.History                         # Returns the MesaLogs object this Profile is a part of.
Profile.HistoryState('star_age')        # Returns the 'star_age' value at this profiles' model number in the MesaLogs object
Profile.HistoryState('all')             # Returns the entire row for this profile's model number in the MesaLogs object
```