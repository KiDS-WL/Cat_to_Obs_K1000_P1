###############
# @file ldac.py
# @author Douglas Applegate & Thomas Erben
# @date 27/11/2013
#
# @brief Utilities to make accessing LDAC cats easier within Python
###############

# HISTORY INFORMATION:
# ====================
#
# 01.09.2010:
# I included treatment of slices through vectorkeys in
# LDACCat.__getitem__
#
# 09.09.2010:
# I made the module more robust against the non-existence
# of necessary libraries
#
# 28.09.2010:
# Bug fix: The from ... import statement must appear on the
# top of the file.
#
# 20.11.2013:
# significant extensions to the code; first proper implementations
# of LDAC tables and LDAC catalogue objects.
#
# 25.11.2013:
# - In the LDACCat class the header of the original catalogue is
#   preserved when storing the catalogue to a file.
# - I add a method to add HISTORY keywords to the image header of
#   the catalogue
#
# 27.11.2013:
# In the LDACTable I replaced the test for an empty LDAC table. Now
# it is done by testing the data element for 'None'. Before it was
# done by the 'size' method on the HDU. It seems that this size element
# of HDU changed from a method to a simple 'int' in different pyfits
# versions. Hence it is useless for code that should be compatible
# to different versions of pyfits.
#
# 28.02.2014:
# In the LDACTable 'setitem' method we do no longer create a new table
# immediatly after a new column was added. This new table creation
# takes long for big catalogues and it leads to necessary operations
# if we add many columns without doing any operations in between. We
# therefore only update column definitions and raise a flag that this
# table needs an 'update'. The update then is done whenever necessary.
#
# 28.03.2014:
# - I corrected function '__setitem__' in LDACCat. The previous
#   implementation was untested and did not work at all!
# - I transfered Dominiks implementation of of the 'add' function in the
#   LDACTable. It enables concatenation of two tables with the same
#   keys and key types.
#
# 24.07.2015:
# - The pyfits module was replaced by the equivalent 'astrpy.io.fits'.
#   pyfits will not be supported anymore at some point in the future.
# - We substituted the depreceated astropy function 'new_table' with
#   the new BinTableHDU.from_columns. The former will no longer
#   be supported in a future version of astropy.
# - We substituted constructs such as 'self.hdu.data == None' to
#   'self.hdu.data is None' to avoid a 'FutureWarning: comparison to
#   `None` will result in an elementwise object comparison in the future.'
#   warning meesage. Although the expected future change of the '=='
#   operator would have no effect in the special cases the warning message
#   is confusing!
# - We do no longer provide detailed error messages if the numpy or astropy
#   (former pyfits) modules are missing. They are quite standard and
#   sufficiently known by now.
#
# 16.08.2015:
# - Bug fix: The LDACcat saveas function did not check whether tables
#   in the catalogue not to be updated before writing the whole catalogue
#   to a file - fixed!
# - In the LDACTable class I had to change the name of the private
#   function '__update()' to '_update()' to be able to use it within
#   the LDACCat class.
#
# 09.11.2015:
# I implemented the __delitem__ method in the LDACTable class
#
# 28.07.2016:
# I converted the script to python3.
#
# 28.06.2017:
# The __add__ function was moved from the LDACCat to the LDACTable object.
# I never used this function and it was added by Dominik. I wrongly
# transfered his code to the wrong class at the time.
#
# 18.05.2018:
# substituted clobber argument with overwrite (astropy deprecation
# warning)
#
# 18.06.2018:
# I modified the add_history function to add a date/time stamp when the history
# was added.
#
# 03.09.2018:
# - I added handling of comments from keys:
#   There currently is not standardised way to handle comments for table keys
#   and astropy deleted those from non-standard conventions if new tables are
#   craeted. We now 'manually' create headers with comments in the TCOMM
#   notation and as comments from the TTYPE keyword.
# - I added primitive treatment of units (you can set and retrieve them for
#   keys but not yet do calculations with them)
#
# 05.09.2018:
# Bug fix: The change from 03.09. introduced a bug which did not allow
# creation of new tables.

"""
Wrapper module to work with LDAC catalogues and tables
"""



# standard-library includes:
import sys
import datetime as dt
import astropy.io.fits as aif
import numpy as np

class LDACCat(object):
    """
    Class to represent an LDAC catalogue
    """

    def __init__(self, cat=None):
        """
        An LDAC catalogue can be instantiated either as an empty catalogue
        or with an existing catalogue on disk.

        >>> a = ldac.LDACCat('mag.cat') # reads the catalogue 'mag.cat' into
                                        # the variable 'a'.
        """

        # The LDACCcat object contains a list of LDAC tables.  We
        # internally also keep the header of the PrimaryHDU. It is
        # reused when the catalogue is saved to a file.

        # for an empty catalogue this list is empty:
        self.ldactables = []
        self.header = None

        if cat != None:
            # read tables from a catalogue on disk:
            if type(cat) == type("a"):
                hdulist = aif.open(cat)

                for hdu in hdulist:
                    if isinstance(hdu, aif.PrimaryHDU) == True:
                        self.header = hdu.header
                    if isinstance(hdu, aif.BinTableHDU) == True:
                        self.ldactables.append(LDACTable(hdu))

    def __len__(self):
        """
        return the number of LDAC tables in this catalogue

        >>> b = len(a)  # number of LDAC tables in catalogue 'a'.
        """

        return len(self.ldactables)

    def __getitem__(self, tablename):
        """
        returns the named LDAC table. Returns 'None' if the table does
        not exist.

        Example:
        >>> b = a['OBJECTS'] # returns in 'b' the LDAC table with name
                             # 'OBJECTS' from catalogue 'a'
        """

        result = None
        for table in self.ldactables:
            if table.hdu.name == tablename:
                result = table

        return result

    def __setitem__(self, name, table):
        """
        adds or replaces an LDAC table in this catalogue

        >>> a['NEW_TABLE'] = b['OBJECTS'] # adds the new table 'NEW_TABLE' in
                                          # 'a' from table 'OBJECTS' in 'b'.
        """

        if isinstance(table, LDACTable):
            # check whether a table with name exists already:
            exists = False

            for i in range(len(self.ldactables)):
                if self.ldactables[i].hdu.name == name:
                    self.ldactables[i] = table
                    exists = True

            if exists == False:
                table.setname(name)
                self.ldactables.append(table)

    def tables(self):
        """
        returns the names of the contained LDAC tables

        >>> c = a.tables()  # gives a list of table names in catalogue 'a'
        """
        tablenames = []

        for table in self.ldactables:
            tablenames.append(table.hdu.name)

        return tablenames

    def __iter__(self):
        return self.ldactables.__iter__()

    def __contains__(self, tablename):
        """
        check whether a table with name 'tablename' is present
        in this catalogue
        """

        return tablename in self.tables()

    def has_table(self, tablename):
        """
        check whether a table with name 'tablename' is present
        in this catalogue

        >>> c = a.has_table('OBJECTS') # returns 'True' if a table named
                                       # 'OBJECTS' is in catalogue 'a'
        """

        return self.__contains__(tablename)

    def add_history(self, history):
        """
        add a history keyword to the header of the catalogue.
        In addition to the history string itself, a date/time
        stamp is added.

        >>> a.add_history('Catalogue created')
        """

        # create an empty header if necessary
        if self.header is None:
            self.header = aif.Header()

        curr_date = dt.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

        # just delegate the work to an astropy method:
        self.header.add_history('') # empty line for separation from other
                                    # comment/history lines
        entry = "history entry added at: %s" % (curr_date)
        self.header.add_history(entry)
        self.header.add_history(history)

    def saveas(self, filename, checksum=False, overwrite=False):
        """
        save the LDAC catalogue to a file.

        if overwrite=True an existing file is overwritten.

        >>> a.saveas('test.cat') # saves LDAC catalogue 'a' with all its
                                 # tables to file 'test.cat'
        """

        primaryHDU = aif.PrimaryHDU(header=self.header)
        hdulist = aif.HDUList([primaryHDU])

        for table in self.ldactables:
            if table.update == 1:
                table._update()
            hdulist.append(table.hdu)

        hdulist.writeto(filename, checksum=checksum, overwrite=overwrite)


class LDACTable(object):
    """
    Class to represent an LDAC table
    """

    def __init__(self, hdu=None):
        """
        An LDAC table can be instantiated either as am empty table
        or with an astropy BinaryTable HDU (existing table).
        """

        # dictionaries to keep comments and units
        self._key_comments = {}
        self._key_units = {}
        self._key_ucd = {}
        self._key_tindx = {}

        if hdu is None:
            self.hdu = aif.BinTableHDU()
            self.hdu.data = None

            # We make sure that the table has 'some' proper name:
            self.hdu.name = "DEFAULT"
        else:
            self.hdu = hdu

            # collect comments and units from all the table keys:
            for i in range(len(hdu.columns.names)):
                name = hdu.header.get('TTYPE%d' % (i + 1))

                # first try comment from TCOMM and then comment from TTYPE:
                try:
                    curr_comm = hdu.header['TCOMM%d' % (i + 1)]
                except:
                    curr_comm = hdu.header.comments('TTYPE%d' % (i + 1))
                    
                self._key_comments[name] = curr_comm
                
                # retrieve unit if available
                try:
                    self._key_units[name] = hdu.header['TUNIT%d' % (i + 1)]

                except:
                    pass
                                            
        self.update = 0 # does the table need an update (e.g. when
                        # new columns were added?)

    def __len__(self):
        """
        return the number of table entries (objects)
        """

        if self.update == 1:
            self._update()

        # 'self.hdu.data' leads to an exception for an empty catalogue.
        # Hence we check for this first:
        if self.hdu.size == 0:
            return 0
        else:
            return len(self.hdu.data)

    def __getitem__(self, key):
        """
        returns the contents of an existing LDAC key as numpy array

        Example:
        >>> b = a['Xpos'] # store in 'b' the contents (numpy array)
                          # of key 'Xpos' from table 'a'.
        """

        if self.update == 1:
            self._update()

        if type(key) == type(5) or \
                type(key) == type(slice(5)):
            return self.hdu.data[key]

        if type(key) == type("a"):
            # we need to deal with slices through vector keys
            # such as 'MAG_APER(2)'
            startind = key.find("(")
            endind = key.find(")")

            if startind > 0 and endind > 0:
                keyname = key[:startind]
                keyindex = int(key[startind + 1:endind]) - 1

                try:
                    return self.hdu.data.field(keyname)[:, keyindex]
                except AttributeError:
                    raise KeyError(key)
            else:
                try:
                    return self.hdu.data.field(key)
                except AttributeError:
                    raise KeyError(key)

        raise TypeError

    def __setitem__(self, key, val):
        """
        set values of an LDAC table

        a['Xpos'] = b # sets the key 'Xpos' in the table 'a' to the
                      # values in numpy array 'b'. If the key does
                      # not yet exist, it is created.
        """
        # VERY uncomplete implementation for the moment!
        # - we only treat scalars for the moment!
        # - we do not check whether the key types match
        #   when an existing key is overwritten

        # sanity checks: the column name must be a string and
        # the value arrays length must match the table data
        # dimension:
        if type(key) == type("a"):
            # The first condition applies to an empty table:
            if self.hdu.data is None or len(val) == self.hdu.data.size:
                # If necessary add a new column to the table
                if self.__contains__(key) == True:
                    # quite some things might go wrong here
                    # (same data type, etc.)
                    self.hdu.data.field(key)[:] = val
                else:
                    # determine format for the new column:
                    colformat = ""
                    if np.issubdtype(val.dtype, float) == True:
                        colformat = "1E"

                    if np.issubdtype(val.dtype, int) == True:
                        colformat = "1I"

                    if np.issubdtype(val.dtype, np.string_) == True or \
                       np.issubdtype(val.dtype, np.unicode_):
                        string_length = val.itemsize
                        colformat = "%dA" % (string_length)

                    # now create the new column and create a 'new' table
                    # with the old plus the new column (I did not find a
                    # way to just append a new column to an existing
                    # table!):
                    newcolumn = aif.Column(name=key, format=colformat,
                                              array=val)

                    #  If you want the new columns appended at the top of the file (e.g an ID) then do it like this
                    self.hdu.columns = aif.ColDefs([newcolumn]) + self.hdu.columns
                    #  If you want the new columns appended at the bottom of the file then do it like this
                    #self.hdu.columns += aif.ColDefs([newcolumn]) 
                    self._key_comments[key] = ""
                    self.update = 1

        #raise NotImplementedError

    def __delitem__(self, key):
        if self.__contains__(key):
            self.hdu.columns.del_col(key)
            del self._key_comments[key]
            del self._key_units[key]
            self.update = 1


    def __add__(self, b):
        """
        Appends table b to table a and returns a new LDAC table.
        Tables 'a' and 'b' must be identical keys and key types.

        >>> c = a + b   # appends table b to a and saves it
                        # as a LDAC table again
        """
        # First check if both tables have the same number of
        # columns:
        if len(list(self.keys())) != len(list(b.keys())):
            print("Tables do not have the same number of columns / keywords!")
            print("First table has " + str(len(list(self.keys()))) + \
                " colums / keywords.")
            print("Second table has " + str(len(list(b.keys()))) + \
                " colums / keywords.")
            return None

        # Now let's check if all keywords from the first table are also
        # present in the second table and also at the same place!
        for i in range(len(list(self.keys()))):
            if (list(self.keys())[i] in b) == False:
                print("Key " + str(list(self.keys())[i]) + \
                    " is not present in the second table!")
                return None

        selfrows = self.hdu.data.shape[0]
        brows = b.hdu.data.shape[0]
        nrows = selfrows + brows
        hdu = aif.BinTableHDU.from_columns(self.hdu.columns, nrows=nrows)
        hdu = self-__correct_header(hdu)

        for i in list(self.keys()):
            hdu.data.field(i)[:selfrows] = self.hdu.data.field(i)
            hdu.data.field(i)[selfrows:] = b.hdu.data.field(i)

        hdu.header = self.hdu.header
        hdu.header['NAXIS2'] = nrows
        hdu.columns = self.hdu.columns
        hdu.name = self.hdu.name

        return LDACTable(hdu)

    def __correct_header(self, hdu):
        """
        'corrects' the header for a newly created binary table.
        key comments are not standardized and hence astropy deleted
        'non-standardized' comments. We add here comments with the
        TCOMM convention as as comments of the TTYPE keyword.

        input:
        - hdu: the HDU of a newly created binary FITS table from
               BinTableHDU
        return:
          The original HDU with the 'corrected' header.
        """

        for i in range(len(hdu.columns.names)):
            name = hdu.header.get('TTYPE%d' % (i + 1))

            if name in self._key_comments and \
               len(self._key_comments[name]) > 0:
                hdu.header.comments['TTYPE%d' % (i + 1)] = \
                  self._key_comments[name]
                hdu.header.set('TCOMM%d' % (i + 1), self._key_comments[name],
                               after='TTYPE%d' % (i + 1))

            if name in self._key_units and \
               len(self._key_units[name]) > 0: 
                hdu.header.set('TUNIT%d' % (i + 1), self._key_units[name],
                               after='TTYPE%d' % (i + 1))
                               
            if name in self._key_ucd and \
               len(self._key_ucd[name]) > 0:
                hdu.header.set('TUCD%d' % (i + 1), self._key_ucd[name],
                               after='TTYPE%d' % (i + 1))

            if name in self._key_tindx and \
               len(self._key_tindx[name]) > 0:
                hdu.header.set('TINDX%d' % (i + 1), self._key_tindx[name],
                               after='TTYPE%d' % (i + 1))

        
        return hdu

    def _update(self):
        # update the table if necessary:
        newtabhdu = aif.BinTableHDU.from_columns(self.hdu.columns)
        newtabhdu = self.__correct_header(newtabhdu)
        newtabhdu.name = self.hdu.name
        self.hdu = newtabhdu
        self.update = 0

    def get_comment(self, key):
        """
        return the comment of a key

        >>> tab.get_comment('x')  # returns comment from key 'x' in table
                                  # tab as a string
        """

        if key in self._key_comments[key]:
            return self._key_comments[key]
        else:
            return ""

    def set_comment(self, key, comm):
        """
        set comment of a key

        >>> tab.set_comment('x', 'x position')

        """

        self._key_comments[key] = comm

    def set_unit(self, key, unit):
        """
        set the unit of a key

        >>> tab.set_unit('x', 'pix')
        """

        self._key_units[key] = unit
        self.hdu.columns[key].unit = unit
        
    def set_ucd(self, key, ucd):
        """
        set the UCD of a key

        >>> tab.set_UCD('x', 'meta.code')
        """

        self._key_ucd[key] = ucd

    def set_tindx(self, key, tindx):
        """
        set the TINDX of a key

        >>> tab.set_TINDX('x', 'T')
        """

        self._key_tindx[key] = tindx


    def keys(self):
        """
        returns the names of the keys contained in this table

        >>> b = a.keys() # store a list of keynames of table 'a' in
                         # 'b'.
        """

        if self.update == 1:
            self._update()

        return self.hdu.columns.names

    def __iter__(self):
        if self.update == 1:
            self._update()

        return self.hdu.data.__iter__()

    def __contains__(self, item):
        if self.update == 1:
            self._update()

        return item in list(self.keys())

    def has_key(self, key):
        """
        tests whether the table contains a certain key.

        >>> b = a.haskey('Xpos') # returns 'True' if table 'a' contains
                                 # a key with name 'Xpos'.
        """

        return self.__contains__(key)

    def filter(self, mask):
        if self.update == 1:
            self._update()

        newtable = aif.BinTableHDU(data=self.hdu.data[mask],
                                   header=self.hdu.header)
        newtable = self.__correct_header(newtable)

        return LDACTable(newtable)

    def setname(self, name):
        """
        set/change the name of the LDAC table.

        >>> a.setname('TESTTABLE') # set/change the name of the LDAC table
                                   # in 'a' to 'TESTTABLE'.
        """

        self.hdu.name = name

    def saveas(self, filename, checksum=False, overwrite=False):
        """
        save the LDAC table as a catalogue. The produced
        catalogue will only consist of this table!

        overwrite=True overwrites an existing file with the
        new catalogue

        >>> a.saveas('table.cat') # saves the LDAC table in 'a'
                                  # to file 'table.cat'
        """

        if self.update == 1:
            self._update()

        self.hdu.writeto(filename, checksum=checksum, overwrite=overwrite)


def openObjects(hdulist, table='OBJECTS'):
    tablehdu = None
    for hdu in hdulist:
        # In a regular LDAC catalogue the primary header
        # does not have an EXTNAME keyword and 'hdu.header['EXTNAME']'
        # leads to a KeyError exception which we just ignore:
        try:
            if table == hdu.header['EXTNAME']:
                tablehdu = hdu
        except KeyError:
            pass

    if tablehdu is None:
        print("Table %s not present in catalogue %s" % (table,
                                                        hdulist.filename()))
        print("Creating an empty LDAC table")

    return LDACTable(tablehdu)

def openObjectFile(filename, table='OBJECTS'):
    hdulist = aif.open(filename)
    if hdulist is None:
        return None

    return openObjects(hdulist, table)
