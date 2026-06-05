'''McMesa.py

A wrapper for mesa_reader that provides alternate functionality to the standard mesa_reader library.
Loads a mesa history file and automatically parses & links models to profiles.
Profiles are loaded as numpy ndarrays and have custom defined methods to make accessing slightly simpler.

'''

import mesa_reader as mr
import numpy as np

class ProfileWrapper:
    profile:mr.MesaData
    ModelNumber:int
    ProfileNumber:int
    Priority:int
    History:MesaLogs
    
    _path:str
    _loaded:bool = False

    ColumnNames:tuple[str] # just has list of columns
    Columns:dict # stores col_name:idx pairings
    Data:np.ndarray # has actual data
    Shape:tuple # tuple

    AutoLower:bool = False
    AutoZone:bool = False
    AutoLog:bool = True

    def __init__(self, root:MesaLogs, path:str, mnum:int, pnum:int, priority:int, load:bool = True):
        self.ModelNumber = mnum
        self.ProfileNumber = pnum
        self.Priority = priority
        self._path = path
        self.History = root

        if load:
            self._load()
        # print(self.Data[:,0]) access ROW 0 (zone)
        # print(self.Data[0,:]) # access COLUMN 0 (data)

    def _load(self):
        if self._loaded:
            return False
        self.profile = mr.MesaData(self._path)
        # self.ColumnNames = self.profile.bulk_names
        # self.Columns = {}
        # tmp_data = []
        # n = 0
        # for col in self.ColumnNames:
        #     self.Columns[col] = n
        #     n += 1
        #     tmp_data.append(self.profile.data(col))

        # self.Data = np.array(tmp_data)
        # self.Shape = self.Data.shape

        self.Data, self.Shape, self.ColumnNames, self.Columns = _Conv_Mesa_NDArray(self.profile)
        self._loaded = True
        return True

    def __getitem__(self, key):
        return _GetItem(self.Data, key, self.ColumnNames, self.Columns, self.AutoZone, self.AutoLog, self.AutoLower, lambda k : k - 1)
        
    def __setitem__(self, key, value):
        if type(key) is tuple:
            a, b = key #data, zone

            if self.AutoZone:
                b = b-1

            if type(a) is str:
                # TODO log check
                
                if self.AutoLower:
                    a = a.lower()
                # print(f'{a}, {b} = {value}')
                self.Data[self.Columns[a],b] = value
            else:
                self.Data[a, b] = value
        else:
            _accepted_types = [list, np.ndarray]

            def check_arr(obj, ndarr):
                if ndarr.shape[0] != obj.shape[1]:
                    raise ValueError
                return ndarr
            
            _fix_types = {
                list : lambda x : np.array(x),
                np.ndarray : lambda x : check_arr(self.Data, x)
            }
            if type(value) not in _accepted_types:
                raise TypeError(value)
            
            value = _fix_types[type(value)](value)
            
            if type(key) is str:
                # TODO log check
                if self.AutoLower:
                    key = key.lower()
                self.Data[self.Columns[key]] = value
            else:
                self.Data[key] = value
    
    def AddColumn(self, column:str, contents=None, allowOverwrite:bool = False):
        '''
        Prepends a column to the profile data. If the contents are not specified, defaults to zeroes.

        Contents can be either a numpy ndarray that has the shape (1, Data.shape[1]) (i.e. 1 axis, as long as any column in this profile)
        or a function that takes a single parameter for a given zone dict and returns a value for the contents of that column.
        '''

        if column in self.ColumnNames and not allowOverwrite: # the column already exists
            raise KeyError(column)
        
            # default = np.zeros(self.Data.shape)
            # tmp_data = np.append(default, self.Data, 1)
            # self.Data = self.profile.bulk_data = tmp_data
            # self.Columns = self.profile.bulk_names = [column] + self.profile.bulk_names
        if contents is not None:
            # catches if its a function
            if callable(contents):
                max_idx = int(self['zone'][-1]+1)
                temp = np.array([contents(self.Zone(e)) for e in range(1, max_idx)])
                contents = [temp]
            elif contents.shape != (1, self.Data.shape[1]):
                raise ValueError(contents)
        else:
            contents = np.zeros((1, self.Data.shape[1]), dtype=np.dtypes.Float64DType)


        if self.AutoLower:
            column = column.lower()

        self.Data = np.append(self.Data, contents, axis=0)
        self.ColumnNames = np.append(self.ColumnNames, column)
        self.Shape = self.Data.shape
        self.Columns[column] = len(self.ColumnNames)-1
        return

    def Zone(self, zone):
        '''Returns a dict for the data at a particular zone.'''
        return dict(zip(self.ColumnNames, self.Data[:,zone-1]))
    
    def ConstructNewData(self, new_columns:str|list[str], root_column:str, base_function:function, step_function:function, direction='asc'):
        '''
        Constructs new data into one or more new column(s) based on the contents of the rest of the profile. Sorts by a root column in either ascending or descending order, applying a base function lambda at the first row out of the sorted view of root_column, and then a step function at every column after.

        The base function lambda should take one parameter for that particular zone's dict.
        The step function lambda should take two parameters - the current zone dict and the prior zone dict.

        Both lambdas should return as many values as new columns to create. 

        '''
        if type(new_columns) is not list:
            new_columns = [new_columns]

        for col in new_columns:
            self.AddColumn(col)

        root = self[root_column]
        indexes = np.argsort(root)
        if direction == 'desc':
            indexes = np.flip(indexes)
        
        # Apply base_function to 0
        zone_0 = self['all', indexes[0]]
        values = base_function(zone_0)
        if type(values) is not list and type(values) is not tuple:
            values = [values]
        for col, val in zip(new_columns, values):
            self[col, indexes[0]] = val

        # Iterate through indexes from 1 onwards and apply base_function.
        for idx in range(1, len(indexes)):
            cur_zone = self['all', indexes[idx]]
            prev_zone = self['all', indexes[idx-1]]

            values = step_function(cur_zone, prev_zone)
            if type(values) is not list and type(values) is not tuple:
                values = [values]

            for col, val in zip(new_columns, values):
                self[col, indexes[idx]] = val

    def HistoryState(self, key:str):
        return self.History[key, self.ModelNumber]
class MesaLogs:
    import numpy as np

    class ProfileSet:
        class ProfileIterator:
            _idx:int
            _pset:MesaLogs.ProfileSet

            def __init__(self, pset:MesaLogs.ProfileSet):
                self._idx = 0
                self._pset = pset

            def __iter__(self):
                return self
            
            def __next__(self):
                if self._idx < len(self._pset.ModelNumbers):
                    entry = self._pset[self._pset.ModelNumbers[self._idx]]
                    self._idx += 1
                    return entry
                raise StopIteration

        _profiles:dict # hashmap of index-wrapper sets
        ModelNumbers:tuple[int,...]
        LoadCount:int
        Count:int
        
        def __init__(self, root:MesaLogs, profile_index, path, autoLoad:bool = True, selection_profilewraps=None):
            if selection_profilewraps != None:
                self._profiles = {}
                temp_mnums = []
                self.LoadCount = 0
                for selection in selection_profilewraps:
                    mnum = selection.ModelNumber
                    temp_mnums.append(mnum)
                    self._profiles[mnum] = selection
                    if selection._loaded:
                        self.LoadCount += 1
                    self.ModelNumbers = tuple(temp_mnums)
                self.Count = len(selection_profilewraps)
            else:
                self._profiles = {}
                self.ModelNumbers = profile_index.model_numbers
                self.LoadCount = 0
                self.Count = len(profile_index.profile_numbers)
                for pidx in profile_index.profile_numbers:
                    mnum = profile_index.model_with_profile_number(pidx)
                    priority = 0
                    self._profiles[mnum] = ProfileWrapper(root, f'{path}/profile{pidx}.data', mnum, pidx, priority, load=autoLoad)

                if autoLoad:
                    self.LoadCount = len(self.ModelNumbers)

        def __iter__(self):
            return MesaLogs.ProfileSet.ProfileIterator(self)
        
        def __getitem__(self, key):
            if type(key) is slice:
                # Create slice representation of this profile set... yesh
                # figure out all within the range of the slice
                lo, hi, step = key.start, key.stop, key.step
                lo = np.min(self.ModelNumbers) if lo is None else lo
                hi = np.max(self.ModelNumbers)+1 if hi is None else hi
                step = 1 if step is None else step

                if lo >= hi or step <= 0:
                    raise ValueError(f'{lo} < {hi} AND {step} > 0 must be true!')

                # create a list of temp nums 
                temp_mnums = np.array(self.ModelNumbers)
                select_mnums = temp_mnums[(temp_mnums >= lo) & (temp_mnums <= hi)]

                select_profilewraps = []
                for k in range(0, len(select_mnums), step):
                    selection = select_mnums[k]
                    select_profilewraps.append(self._profiles[selection])
                return MesaLogs.ProfileSet(None, None, None, None, select_profilewraps)
            elif key in self._profiles:
                if self._profiles[key]._load():
                    self.LoadCount += 1
                return self._profiles[key]
            raise KeyError(key)
        
    name:str
    path:str

    history:mr.MesaData
    _modelmap:dict

    Data:np.ndarray
    Shape:tuple[int, int]
    ColumnNames:tuple[str,...]
    Columns:dict

    Profiles:ProfileSet
    ModelProfiles:tuple[int,...]

    AutoCorrectModelNum:bool
    
    def __init__(self, path:str, autoLoad:bool = True):
        splitchar = '/' if '/' in path else '\\'
        if path.endswith(splitchar):
            path = path[:-1]
        self.name = path.split(splitchar)[-1][len('LOGS_'):]
        self.path = path

        self.history = mr.MesaData(f'{self.path}/history.data')
        # self.Columns = self.history.bulk_names

        profile_index = mr.MesaProfileIndex(f'{self.path}/profiles.index')
        self.ModelProfiles = profile_index.model_numbers

        self.Profiles = MesaLogs.ProfileSet(self, profile_index, path, autoLoad)

        self.Data, self.Shape, self.ColumnNames, self.Columns = _Conv_Mesa_NDArray(self.history)
        self.AutoCorrectModelNum = True

        model_nums = _GetItem(self.Data, 'model_number', self.ColumnNames, self.Columns, False, False, False, None)
        # create reverse hashmap
        self._modelmap = {}
        for k in range(len(model_nums)):
            self._modelmap[model_nums[k]] = k

        column='hasProfile'
        contents = np.zeros(self.Data.shape[1], dtype=np.dtypes.Float64DType)
        # print(f'{contents.shape}\n{self.Data.shape}')
        for index in self.ModelProfiles:
            # print(f'{index} => {self._modelmap[index]}')
            contents[self._modelmap[index]] = 1
        
        contents = [contents]

        self.Data = np.append(self.Data, contents, axis=0)
        self.ColumnNames = np.append(self.ColumnNames, column)
        self.Shape = self.Data.shape
        self.Columns[column] = len(self.ColumnNames)-1


    def __getitem__(self, key):
        return _GetItem(self.Data, key, self.ColumnNames, self.Columns, self.AutoCorrectModelNum, True, False, self.GetModelNumberIndex)
    
    def GetModelNumberIndex(self, number):
        return self._modelmap[number]
        
def _Conv_Mesa_NDArray(profile)->tuple[np.ndarray, tuple[int, int], list[str], dict]:
    ColumnNames = profile.bulk_names
    Columns = {}
    tmp_data = []
    n = 0
    for col in ColumnNames:
        Columns[col] = n
        n += 1
        tmp_data.append(profile.data(col))

    Data = np.array(tmp_data)
    Shape = Data.shape

    return Data, Shape, ColumnNames, Columns

def _GetItem(Data, key, ColumnNames, Columns, AutoCorrectPrimary, AutoLog, AutoLower, PrimaryCorrector:function):
    if type(key) is tuple:
        a, b = key #data, zone OR model number

        if AutoCorrectPrimary:
            b = PrimaryCorrector(b)

        if a is None or a.lower() == 'all':
            return dict(zip(ColumnNames, Data[:,b]))
        

        if type(a) is str:
            # TODO log check
            if AutoLower:
                a = a.lower()
            
            if AutoLog:
                # Try see if there exists either 'log_{a}' or 'log{a}' in columns: if there is, then change a to that match and apply as exponent
                if f'log_{a}' in Columns:
                    return np.pow(10, Data[Columns[f'log_{a}'], b])
                if f'log{a}' in Columns:
                    return np.pow(10, Data[Columns[f'log{a}'], b])
            return Data[Columns[a],b]
        else:
            return Data[a, b]
    else:
        if type(key) is str:
            # TODO log check
            if AutoLower:
                key = key.lower()
            if AutoLog:
                # Try see if there exists either 'log_{a}' or 'log{a}' in columns: if there is, then change a to that match and apply as exponent
                if f'log_{key}' in Columns:
                    return np.pow(10, Data[Columns[f'log_{key}']])
                if f'log{key}' in Columns:
                    return np.pow(10, Data[Columns[f'log{key}']])
            return Data[Columns[key],:]
        else:
            return Data[key,:]