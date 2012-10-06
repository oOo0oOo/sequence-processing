import wx
import wx.grid
import os, sys
import cell
import parse
import time

from parse import SequenceError, SequenceWarning

class PreferencesPopup(wx.Dialog):
    '''This is a dialog to change the preferences of a cell, it takes a 
    cell as inpput and parses uses the current parameters.'''
     
    def __init__(self, cell):
        super(PreferencesPopup, self).__init__(None)
        
        pnl = wx.Panel(self)
        sample_modes = ['Standard', 'Normal Distribution']
        self.sample_mode = wx.ComboBox(pnl, choices=sample_modes, 
            style=wx.CB_READONLY)
        
        self.st = wx.StaticText(pnl, label='\nE.g. 0.5 => 0 (50% Probability) or 1 (50%); 3.2 => 3 (80%) or 4 (20%)', size=(350, 65))
        
        q = 'Algorithm to use to sample energy consumption, unit max-action parameter and\n'
        q += 'constitutive promoters translation number (applied every round).\n'
        
        self.sample_label = wx.StaticText(pnl, label=q)
        self.mode_label = wx.StaticText(pnl, label='Sample Mode:   ')
        self.sigma_label = wx.StaticText(pnl, label='    Sigma Factor:   ')
        self.sigma_factor = wx.TextCtrl(pnl, wx.ID_ANY)
        self.sample_mode.Bind(wx.EVT_COMBOBOX, self.OnSelect)
        
        sample_sizer = wx.BoxSizer(wx.HORIZONTAL)
        sample_sizer.Add(self.mode_label)
        sample_sizer.Add(self.sample_mode)
        sample_sizer.Add(self.sigma_label)
        sample_sizer.Add(self.sigma_factor)
        
        self.interact_chk = wx.CheckBox(pnl, -1, "Allow transport into subcells")
        
        self.subs_const = wx.ListCtrl(pnl, -1, size=(200, 100), style=wx.LC_REPORT | wx.SUNKEN_BORDER)
        self.subs_const.InsertColumn(0, 'Substance')
        self.subs_const.InsertColumn(1, 'New per round')
        
        #fill values in
        for sub in cell.subs_const.keys():
            index = self.subs_const.InsertStringItem(sys.maxint, sub)
            self.subs_const.SetStringItem(index, 1, str(cell.subs_const[sub]))
        
        self.subs_const.SetColumnWidth(0, wx.LIST_AUTOSIZE_USEHEADER)
        self.subs_const.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        
        self.add_sub_btn = wx.Button(pnl, wx.ID_ANY, "Add Substance")
        self.remove_sub_btn = wx.Button(pnl, wx.ID_ANY, "Remove Substance")

        self.Bind(wx.EVT_BUTTON, self.add_sub, self.add_sub_btn)
        self.Bind(wx.EVT_BUTTON, self.remove_sub, self.remove_sub_btn)
        
        substance_sizer = wx.BoxSizer(wx.VERTICAL)
        subs_btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        subs_btn_sizer.Add(self.add_sub_btn)
        subs_btn_sizer.Add(self.remove_sub_btn)
        substance_sizer.Add(self.subs_const)
        substance_sizer.Add(subs_btn_sizer)
        
        okbtn = wx.Button(pnl, wx.ID_OK, label='Save Preferences')
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.sample_label)
        main_sizer.Add(sample_sizer)
        main_sizer.Add(self.st)
        main_sizer.Add(self.interact_chk)
        main_sizer.Add(wx.StaticText(pnl, label='\n\nSubstances constantly introduced into top cell:'))
        main_sizer.Add(substance_sizer)
        main_sizer.Add(wx.StaticText(pnl, label='\n'))
        main_sizer.Add(okbtn)
        pnl.SetSizer(main_sizer)
        
        self.SetSize((500, 450))
        self.SetTitle('Preferences: ' + cell.name + ', ID: ' + str(cell.id))
        #self.Centre()
        self.Show(True)
        
        self.current_sigma_factor = cell.normal_factor
        
        if cell.sample_mode == 'standard':
            sample_mode = 'Standard'
        elif cell.sample_mode == 'normal':
            sample_mode = 'Normal Distribution'
        
        self.sample_mode.SetValue(sample_mode)
        self.sigma_factor.SetValue(str(cell.normal_factor))
        
        if cell.interact_with_subcells:
            self.interact_chk.Set3StateValue(wx.CHK_CHECKED)
        
    def add_sub(self, e):
        dlg = wx.TextEntryDialog(None, 'Substance tag (only lowercase characters\nand underscore allowed.','Substance Tag')
        if dlg.ShowModal() == wx.ID_OK:
            tag = dlg.GetValue()
            dlg = wx.TextEntryDialog(None, 'Amount of Molecules added per round','Substance Amount')
            if dlg.ShowModal() == wx.ID_OK:
                amount = dlg.GetValue()
                index = self.subs_const.InsertStringItem(sys.maxint, tag)
                self.subs_const.SetStringItem(index, 1, amount)
    
    def remove_sub(self, e):
        selected = self.subs_const.GetNextSelected(-1)
        if selected != -1:
            self.subs_const.DeleteItem(selected)
    
    def get_preferences(self):
        #Get all entered values
        sample_mode = self.sample_mode.GetValue()
        sigma_factor = float(self.sigma_factor.GetValue())
        
        if sample_mode == 'Standard':
            sample_mode = 'standard'
        elif sample_mode == 'Normal Distribution':
            sample_mode = 'normal'
            
        count = self.subs_const.GetItemCount()
        subs_const = {}
        for row in range(count):
            tag = self.subs_const.GetItem(itemId=row, col=0).GetText()
            number = int(self.subs_const.GetItem(itemId=row, col=1).GetText())
            subs_const[tag] = number
            
        chk_state = self.interact_chk.Get3StateValue()
        if chk_state == wx.CHK_CHECKED:
            interact = 'True'
        elif chk_state == wx.CHK_UNCHECKED:
            interact = 'False'
        
        return sample_mode, sigma_factor, subs_const, interact
        
    def OnSelect(self, e):
        i = e.GetString()
        sig = self.current_sigma_factor
        if i == 'Standard':
            label = '\nE.g. 0.5 => 0 (50% Probability) or 1 (50%); 3.2 => 3 (80%) or 4 (20%)'
        elif i == 'Normal Distribution':
            label = '\nMean = Parameter; Sigma = Sigma_factor * Parameter'
            label += '\n0 < Sigma_factor < 1: Normal use; Sigma_factor > 1: Highly Random'
            if sig == 0:
                sig = 0.2
            
        self.st.SetLabel(label)
        self.sigma_factor.SetValue(str(sig))

class ProgressPopup(wx.Frame):
    def __init__(self, parent, num_steps):
        
        wx.Frame.__init__(self, parent, id=wx.ID_ANY, size = (250, 125), 
                          style=wx.CAPTION)
        self.cur_round = 0
        self.num_steps = num_steps
        self.moving_avg = []
        
        self.gauge = wx.Gauge(self, range=self.num_steps, size=(250, 25))
        self.label = wx.StaticText(self)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.gauge)
        sizer.Add(self.label)
        self.SetSizer(sizer)
        self.Show()
        
        self.gauge.SetValue(0)
        self.before = time.time()
        self.label.SetLabel('\n  Modeling round: {0} / {1} \n'.format(0, self.num_steps))
        
    def end_round(self):
        '''Update the gauge and time estimate'''
        
        self.cur_round += 1
        
        time_round = round(time.time() - self.before,4)
        self.before = time.time()
        
        self.moving_avg.append(time_round)
        if len(self.moving_avg) > 8:
            self.moving_avg.pop(0)
        
        cur_avg = sum(self.moving_avg)/len(self.moving_avg)
        if cur_avg > 0.1:
            cur_avg = time_round
        
        cur_estimate = cur_avg * (self.num_steps - self.cur_round)
                
        if cur_estimate <= 1:
            cur_estimate = '< 1'
        
        elif 1 < cur_estimate <= 5:
            cur_estimate = '< 5'
            
        elif 5 < cur_estimate <= 10:
            cur_estimate = '< 10'
        
        elif 10 < cur_estimate:
            cur_estimate = round(cur_estimate)
            
        message = '\n  Modeling round: {0} / {1} \n  Time needed  for last round: {2}s'.format(self.cur_round + 1, self.num_steps, time_round)
        if self.cur_round != 0:
            message += '\n  Estimated finish at equal speed in {0} s'.format(cur_estimate)
        
        self.label.SetLabel(message)
        self.gauge.SetValue(self.cur_round)

class SequenceTab(wx.Panel):
    """
    Sequence Tab.
    """
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        msizer = wx.BoxSizer(wx.VERTICAL)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        
        box = wx.StaticBox(self, -1, "All Sequences in Editor")
        bsizer = wx.StaticBoxSizer(box, wx.HORIZONTAL)
        
        self.grid = wx.grid.Grid(self, size=(725, 450))
        self.grid.CreateGrid(0, 4)
        
        #Make grid readonly
        grid_attr = wx.grid.GridCellAttr()
        grid_attr.SetReadOnly(True)
        for ind in (2, 3):
            self.grid.SetColAttr(ind, grid_attr)
        
        columns = ['Name', 'Sequence', 'Length', 'Syntax']
        for ind, name in enumerate(columns):
            self.grid.SetColLabelValue(ind, name)
        
        
        bsizer.Add(self.grid)
        
        hsizer.Add((5,10),0,0,0)
        hsizer.Add(bsizer)
        
        msizer.Add((10,10),0,0,0)
        msizer.Add(hsizer)
        self.SetSizer(msizer)
        
        
        self.grid.AutoSizeColumns()
        
    def corners_to_cells(self, top_lefts, bottom_rights):
        """
        From: http://ginstrom.com/scribbles/2008/09/07/getting-the-selected-cells-from-a-wxpython-grid/
        """
    
        cells = []
        for top_left, bottom_right in zip(top_lefts, bottom_rights):
    
            rows_start = top_left[0]
            rows_end = bottom_right[0]
    
            cols_start = top_left[1]
            cols_end = bottom_right[1]
    
            rows = range(rows_start, rows_end+1)
            cols = range(cols_start, cols_end+1)
    
            cells.extend([(row, col)
                for row in rows
                for col in cols])
    
        return cells 

    def get_selected_cells(self):
        """
        From: http://ginstrom.com/scribbles/2008/09/07/getting-the-selected-cells-from-a-wxpython-grid/
        """
    
        top_left = self.grid.GetSelectionBlockTopLeft()
    
        if top_left:
            bottom_right = self.grid.GetSelectionBlockBottomRight()
            return self.corners_to_cells(top_left, bottom_right)
    
        selection = self.grid.GetSelectedCells()
    
        if not selection:
            row = self.grid.GetGridCursorRow()
            col = self.grid.GetGridCursorCol()
            return [(row, col)]
    
        return selection
        
class CellTab(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        
        
        msizer = wx.BoxSizer(wx.VERTICAL)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        
        #boxes to hold options
        box1 = wx.StaticBox(self, -1, "Current Population")
        b1sizer = wx.StaticBoxSizer(box1, wx.VERTICAL)
        box2 = wx.StaticBox(self, -1, "Selected Cell Details")
        b2sizer = wx.StaticBoxSizer(box2, wx.HORIZONTAL)
        
        #Tree to show all cells and nested subcells
        self.tree = wx.TreeCtrl(self, wx.ID_ANY, size = (375, 450))
        
        #Buttons to expand collapse tree
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.expand_btn = wx.Button(self, wx.ID_ANY, "Expand All")
        self.collapse_btn = wx.Button(self, wx.ID_ANY, "Collapse All")
        
        self.Bind(wx.EVT_BUTTON, self.expand, self.expand_btn)
        self.Bind(wx.EVT_BUTTON, self.collapse, self.collapse_btn)
        
        btnsizer.Add(self.expand_btn)
        btnsizer.Add(self.collapse_btn)
        
        #TextCntrl for Cell details (readonly)
        self.cell_details = wx.TextCtrl(self, wx.ID_ANY, 'No cell selected...',
                       size=(350, 450), style=wx.TE_MULTILINE|wx.TE_READONLY)
        
        b1sizer.Add(self.tree, wx.EXPAND)
        b1sizer.Add(btnsizer)
        b2sizer.Add(self.cell_details)
        
        hsizer.Add((5,10),0,0,0)
        hsizer.Add(b1sizer)
        hsizer.Add((10,10),0,0,0)
        hsizer.Add(b2sizer)
        
        msizer.Add((10,10),0,0,0)
        msizer.Add(hsizer)
        self.SetSizer(msizer)
        
        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnSelChanged, self.tree)
        
        self.current_selection = (0,)
        
    def OnSelChanged(self, event):
        item = event.GetItem()
        if item:
            self.current_adress = self.tree.GetItemPyData(item)
            #get the cell by adress
            cell = self.top_cell.get_subcell_by_adress(self.current_adress)
            self.cell_details.SetValue(cell.parse_status())
            
        event.Skip()
        
    def get_selected(self):
        return self.current_adress
    
    def collapse(self, e):
        self.tree.CollapseAll()
    
    def expand(self, e):
        self.tree.ExpandAll()
        
    def update_tree(self, top_cell):
        self.top_cell = top_cell
        self.tree.DeleteAllItems()
        
        isz = (16,16)
        il = wx.ImageList(isz[0], isz[1])
        fldridx     = il.Add(wx.ArtProvider_GetBitmap(wx.ART_FOLDER,      wx.ART_OTHER, isz))
        fldropenidx = il.Add(wx.ArtProvider_GetBitmap(wx.ART_FILE_OPEN,   wx.ART_OTHER, isz))
        fileidx     = il.Add(wx.ArtProvider_GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, isz))

        self.tree.SetImageList(il)
        self.il = il
        
        #helper function to parse cell data into nice string
        def parse_cell(cell, adress):
            if cell.name == '':
                cell_str = 'Unnamed, '
            else:
                cell_str = cell.name + ', '
            
            cell_str += 'id: ' + str(cell.id) + ', adress: ('
            c = 0
            for ad in adress:
                if c == 1:
                    cell_str += ', '
                cell_str += str(ad)
                c = 1
            cell_str += ')'
            return cell_str
        
        self.tree_items = {}
        
        #set top cell (root item in tree)
        self.tree_items[(0,)] = self.tree.AddRoot(parse_cell(self.top_cell, (0,)))
        self.tree.SetPyData(self.tree_items[(0,)], (0,))
        self.tree.SetItemImage(self.tree_items[(0,)], fldridx, wx.TreeItemIcon_Normal)
        self.tree.SetItemImage(self.tree_items[(0,)], fldropenidx, wx.TreeItemIcon_Expanded)
        self.tree.SelectItem(self.tree_items[(0,)])
        
        def add_sub_cells(cell, previous_adress):
            for ind in cell.cells.keys():
                subcell = cell.cells[ind]
                #change adress
                adress = list(previous_adress)
                adress.append(subcell.id)
                adress = tuple(adress)
                
                #add item to tree
                self.tree_items[adress] = self.tree.AppendItem(self.tree_items[previous_adress], parse_cell(subcell, adress))
                self.tree.SetPyData(self.tree_items[adress], adress)
                self.tree.SetItemImage(self.tree_items[adress], fldridx, wx.TreeItemIcon_Normal)
                self.tree.SetItemImage(self.tree_items[adress], fldropenidx, wx.TreeItemIcon_Expanded)
                
                #call function recursively to catch all cells
                add_sub_cells(subcell, adress)
        
        add_sub_cells(self.top_cell, (0,))
        
        self.tree.ExpandAll()

class SeqGUI(wx.Frame):
    
    def __init__(self, *args, **kwargs):
        super(SeqGUI, self).__init__(size = (800, 600), *args, **kwargs) 
        
        #initialize top cell
        self.top_cell = cell.Cell(name='Top Cell')
        
        #initialize sequence list
        self.sequences = []
        
        #Initialize the User interface
        self.InitUI()
        
        #update the cell information
        self.celltab.update_tree(self.top_cell)
        
    def InitUI(self):

        menubar = wx.MenuBar()
        
        appMenu = wx.Menu()
        seqMenu = wx.Menu()
        cellMenu = wx.Menu()
        
        #Define menu items
        ami = wx.MenuItem(appMenu, wx.ID_ANY, '&About')
        qmi = wx.MenuItem(appMenu, wx.ID_EXIT, '&Quit\tCtrl+W')
        
        lmi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Load sequence (*.txt)')
        smi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Save sequences to file')
        crmi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Create cells from sequences')
        wmi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Enter new sequence')
        cmi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Re-check syntax of sequences')
        comi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Combine sequences')
        rmi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Remove sequences')
        dmi = wx.MenuItem(seqMenu, wx.ID_ANY, '&Sequence details (incomplete)')
        
        rami = wx.MenuItem(cellMenu, wx.ID_ANY, '&Run all cells')
        tmi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Translate sequence in cell')
        splitmi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Split cell')
        dismi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Disintegrate cell')
        
        pmi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Cell preferences')
        ecmi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Edit cell sequence')
        
        exmi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Extract cell sequence')
        savemi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Save cell log to csv')
        
        demi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Delete selected cell')
        fmi = wx.MenuItem(cellMenu, wx.ID_ANY, '&Delete all cells')
        
        #Bind handler functions to click events
        self.Bind(wx.EVT_MENU, self.preferences, pmi)
        self.Bind(wx.EVT_MENU, self.flush_all, fmi)
        self.Bind(wx.EVT_MENU, self.delete_cell, demi)
        self.Bind(wx.EVT_MENU, self.save_cell_log, savemi)
        self.Bind(wx.EVT_MENU, self.read_seq, tmi)
        self.Bind(wx.EVT_MENU, self.split_cell, splitmi)
        self.Bind(wx.EVT_MENU, self.disintegrate_cell, dismi)
        self.Bind(wx.EVT_MENU, self.run_cells, rami)
        self.Bind(wx.EVT_MENU, self.create_cell, crmi)
        self.Bind(wx.EVT_MENU, self.sequence_details, dmi)
        self.Bind(wx.EVT_MENU, self.remove_sequence, rmi)
        self.Bind(wx.EVT_MENU, self.combine_sequence, comi)
        self.Bind(wx.EVT_MENU, self.check_sequence, cmi)
        self.Bind(wx.EVT_MENU, self.enter_sequence, wmi)
        self.Bind(wx.EVT_MENU, self.edit_cell_sequence, ecmi)
        self.Bind(wx.EVT_MENU, self.extract_sequence, exmi)
        self.Bind(wx.EVT_MENU, self.save_file, smi)
        self.Bind(wx.EVT_MENU, self.load_file, lmi)
        
        self.Bind(wx.EVT_MENU, self.OnQuit, qmi)
        self.Bind(wx.EVT_MENU, self.OnAbout, ami)
        
        
        appMenu.AppendItem(ami)
        appMenu.AppendItem(qmi)
        
        
        seqMenu.AppendItem(lmi)
        seqMenu.AppendItem(smi)
        seqMenu.AppendItem(crmi)
        seqMenu.AppendSeparator()
        
        seqMenu.AppendItem(wmi)
        seqMenu.AppendItem(cmi)
        seqMenu.AppendSeparator()
        
        seqMenu.AppendItem(comi)
        seqMenu.AppendItem(rmi)
        seqMenu.AppendSeparator()
        
        seqMenu.AppendItem(dmi)
        
        
        cellMenu.AppendItem(rami)
        cellMenu.AppendItem(tmi)
        cellMenu.AppendSeparator()
        
        cellMenu.AppendItem(splitmi)
        cellMenu.AppendItem(dismi)
        cellMenu.AppendSeparator()
        
        cellMenu.AppendItem(pmi)
        cellMenu.AppendItem(ecmi)
        cellMenu.AppendSeparator()
        
        cellMenu.AppendItem(savemi)
        cellMenu.AppendItem(exmi)
        cellMenu.AppendSeparator()
        
        cellMenu.AppendItem(demi)
        cellMenu.AppendItem(fmi)
        
        menubar.Append(appMenu, '&App')
        menubar.Append(seqMenu, '&Sequence')
        menubar.Append(cellMenu, '&Cell')

        self.SetMenuBar(menubar)
        
        status_panel = wx.Panel(self)
        self.note_book = wx.Notebook(status_panel, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT, size = (800, 500)
                             )

        # Create the sequence tab and add it to the notebook
        self.seqtab = SequenceTab(self.note_book)
        self.seqtab.SetBackgroundColour("White")
        self.note_book.AddPage(self.seqtab, "Sequences")
        
        # Create the pool tab
        self.celltab = CellTab(self.note_book)
        self.celltab.SetBackgroundColour("White")
        self.note_book.AddPage(self.celltab, "Cells")
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.note_book, 1, wx.ALL|wx.EXPAND, 5)
        status_panel.SetSizer(sizer)
        
        self.Layout()
        
        self.SetTitle('Sequence Processing 0.1')
        self.Centre()
        self.Show()

    def add_sequence(self, name, sequence):
        '''add a new sequence to the grid in the sequence tab'''
        parser = parse.Parser()
        parser.sequence = sequence
        #check sequence syntax
        try:
            parser.check_syntax()
            #Change Syntax value
            syntax = 'OK'
                    
        except (SequenceError, SequenceWarning), e:
            #Change Syntax value
            if type(e) == SequenceError:
                syntax = 'Error'
                dial = wx.MessageDialog(None, str(e), 'Sequence Error', 
                                                wx.OK | wx.ICON_ERROR)
                dial.ShowModal()
            elif type(e) == SequenceWarning:
                syntax = 'Warning'
                dial = wx.MessageDialog(None, str(e), 'Sequence Warning', 
                                                wx.OK | wx.ICON_WARNING)
                dial.ShowModal()
        
        #append sequence to sequence table
        next_row = self.seqtab.grid.GetNumberRows()
        self.seqtab.grid.InsertRows(pos=next_row, numRows=1)
        self.seqtab.grid.SetRowLabelValue(next_row, 'Sequence %i'%next_row)
            
        self.seqtab.grid.SetCellValue(next_row, 0, name)
        self.seqtab.grid.SetCellValue(next_row, 1, sequence)
        self.seqtab.grid.SetCellValue(next_row, 2, str(len(sequence)))
        self.seqtab.grid.SetCellValue(next_row, 3, syntax)
        
        self.seqtab.grid.AutoSizeColumns()
    
    def load_file(self, e):
        """ Open a file, check syntax, create new cell"""
        dlg = wx.FileDialog(self, 'Choose a .txt file to load', '', '', '*.txt', wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            file_path = os.path.join(dirname, filename)
            
            #Load File, check syntax
            parser = parse.Parser()
            parser.load_file(file_path)
            parser.clean_sequence()
            
            filename = filename.replace('.txt', '')
            
            #Add sequence
            self.add_sequence(filename, parser.sequence)
            
        dlg.Destroy()
       
    def save_file(self, e):
        
        sel = self.seqtab.get_selected_cells()
        if sel != [(-1, -1)]:
            sequence = ''
            names = []
            for entry in sel:
                sequence += self.seqtab.grid.GetCellValue(entry[0], 1)
                names.append(self.seqtab.grid.GetCellValue(entry[0], 0))
                
            dial = wx.MessageDialog(self, "Would you like to add auto-generated comments to the sequence?", style=wx.YES_NO|wx.CANCEL|wx.ICON_QUESTION)
            res = dial.ShowModal()
            if res == wx.ID_YES:
                verbose = True
            elif res == wx.ID_NO:
                verbose = False
            
            dial.Destroy()
            if res != wx.ID_CANCEL:
                dlg = wx.FileDialog(self, 'Choose a filename', '', '', '*.txt', wx.OPEN)
                if dlg.ShowModal() == wx.ID_OK:
                    filename = dlg.GetFilename()
                    dirname = dlg.GetDirectory()
                    filepath = os.path.join(dirname, filename)
                        
                    #Load File, check syntax
                    parser = parse.Parser()
                    parser.sequence = sequence
                    parser.clean_sequence()
                    try:
                        parser.check_syntax()
                        try:
                            parser.save_file(filepath, names, verbose)
                            #Show success message
                            message = 'Successfully created \"{0}\" at \"{1}\"'.format(filename, dirname)
                            dial = wx.MessageDialog(None, message, 'Info', wx.OK)
                            dial.ShowModal()
                
                        except Exception, e:
                            print Exception, e
                           
                    except (SequenceError, SequenceWarning), e:
                        dial = wx.MessageDialog(None, str(e), 'Sequence Error / Warning', 
                                                            wx.OK | wx.ICON_ERROR)
                        dial.ShowModal()
                        try:
                            parser.save_file(filepath, names, verbose)
                            #Show success message
                            message = 'Successfully created \"{0}\" at \"{1}\"'.format(filename, dirname)
                            dial = wx.MessageDialog(None, message, 'Info', wx.OK)
                            dial.ShowModal()
                        except Exception, e:
                            print Exception, e
                        
                dlg.Destroy()
        else:
            dial = wx.MessageDialog(None, 'No sequences selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
    
    def enter_sequence(self, e):
        dial = wx.TextEntryDialog(self, 'New sequence:', 'Enter new sequence')
        if dial.ShowModal() == wx.ID_OK:
            new_sequence = dial.GetValue()
            #clean all whitespace and comments
            parser = parse.Parser()
            parser.sequence = new_sequence
            parser.clean_sequence()
            
            #Add sequence
            self.add_sequence('', parser.sequence)
            
        dial.Destroy()
    
    def check_sequence(self, e):
        sel = self.seqtab.get_selected_cells()
        
        if sel != [(-1, -1)]:
        
            for entry in sel:
                seq = self.seqtab.grid.GetCellValue(entry[0], 1)
                parser = parse.Parser()
                parser.sequence = seq
                    
                try:
                    parser.check_syntax()
                    #Change Syntax value
                    seq = self.seqtab.grid.SetCellValue(entry[0], 3, 'OK')
                        
                except (SequenceError, SequenceWarning), e:
                    #Change Syntax value
                    if type(e) == SequenceError:
                        seq = self.seqtab.grid.SetCellValue(entry[0], 3, 'Error')
                        dial = wx.MessageDialog(None, str(e), 'Sequence Error', 
                                                    wx.OK | wx.ICON_ERROR)
                        dial.ShowModal()
                    elif type(e) == SequenceWarning:
                        seq = self.seqtab.grid.SetCellValue(entry[0], 3, 'Warning')
                        dial = wx.MessageDialog(None, str(e), 'Sequence Warning', 
                                                    wx.OK | wx.ICON_WARNING)
                        dial.ShowModal()
        else:
            dial = wx.MessageDialog(None, 'No sequences selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
                    
    def combine_sequence(self, e):
        sel = self.seqtab.get_selected_cells()
        sequence = ''
        for entry in sel:
            sequence += self.seqtab.grid.GetCellValue(entry[0], 1)
        
        if len(sequence) > 0:    
            dial = wx.TextEntryDialog(self, 'Name of new sequence:', 'Combine sequences')
            if dial.ShowModal() == wx.ID_OK:
                name = dial.GetValue()
                #clean all whitespace and comments
                parser = parse.Parser()
                parser.sequence = sequence
                parser.clean_sequence()
                
                #Add sequence
                self.add_sequence(name, parser.sequence)
                
            dial.Destroy()
        else:
            dial = wx.MessageDialog(None, 'No sequences selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()

    def remove_sequence(self, e):
        sel = self.seqtab.get_selected_cells()
        if sel != [(-1, -1)]:
            sel = [list(seq) for seq in sel]
            for ind in range(len(sel)):
                to_del = sel[ind][0]
                self.seqtab.grid.DeleteRows(to_del)
                for ind2 in range(ind + 1, len(sel)):
                    if sel[ind2][0] > sel[ind][0]:
                        sel[ind2][0] -= 1
                        
            self.seqtab.grid.AutoSizeColumns
        else:
            dial = wx.MessageDialog(None, 'No sequences selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
            
    def sequence_details(self, e):
        sel = self.seqtab.get_selected_cells()
        if sel != [(-1, -1)]:
            for entry in sel:
                sequence = self.seqtab.grid.GetCellValue(entry[0], 1)
                parser = parse.Parser()
                parser.sequence = sequence
                
                border = '###################################################'
                print border + '\n' + border + 'Sequence Details:'
                print parser.sequence_info()  
        else:
            dial = wx.MessageDialog(None, 'No sequences selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
            
    def create_cell(self, e):        
        sel = self.seqtab.get_selected_cells()
        if sel != [(-1, -1)]:
            init = False
            dial = wx.MessageDialog(self, "Would you like to initialize the cells\nby translating sequence from position 0?", style=wx.YES_NO|wx.CANCEL|wx.ICON_QUESTION)
            res = dial.ShowModal()
            if res == wx.ID_YES:
                init = True
            dial.Destroy()
            
            if res != wx.ID_CANCEL:
                for entry in sel:
                    seq = self.seqtab.grid.GetCellValue(entry[0], 1)
                    try:
                        parser = parse.Parser()
                        parser.sequence = seq
                        parser.check_syntax()
                        #if everything is ok create cell
                        name = self.seqtab.grid.GetCellValue(entry[0], 0)
                        top_adress = self.celltab.get_selected()
                        cell = self.top_cell.get_subcell_by_adress(top_adress)
                        new_id = cell.new_sub_cell(sequence = seq, name = name)
                        if init:
                            adress = list(top_adress)
                            adress.append(new_id)
                            self.top_cell.read_seq(start = 0, adress = tuple(adress))
                                
                    except (SequenceError, SequenceWarning), e:
                        #Show sequence error dialog
                        if type(e) == SequenceError:
                            dial = wx.MessageDialog(None, str(e), 'Could not create cell', 
                                                            wx.OK | wx.ICON_ERROR)
                            dial.ShowModal()
                            
                        elif type(e) == SequenceWarning:
                            #also create cell if sequence warning was triggered
                            name = self.seqtab.grid.GetCellValue(entry[0], 0)
                            adress = self.celltab.get_selected()
                            cell = self.top_cell.get_subcell_by_adress(adress)
                            new_id = cell.new_sub_cell(sequence = seq, name = name)
                            if init:
                                adress = list(top_adress)
                                adress.append(new_id)
                                self.top_cell.read_seq(start = 0, adress = tuple(adress))
                            
                            #Show warning dialog
                            dial = wx.MessageDialog(None, str(e), 'Created cell', 
                                                            wx.OK | wx.ICON_WARNING)
                            dial.ShowModal()
        else:
            dial = wx.MessageDialog(None, 'No sequences selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
                  
        self.celltab.update_tree(self.top_cell)
            
    def run_cells(self, e):
        dial = wx.TextEntryDialog(self, 'Number of rounds to run the cells:', 'Run Cells')
        dial.SetValue("1")
        if dial.ShowModal() == wx.ID_OK:
            num_rounds = int(dial.GetValue())
            
            #start up progress bar
            progress = ProgressPopup(None, num_rounds)
                
            for i in range(num_rounds):
                self.top_cell.run_all()
                progress.end_round()
                
            self.celltab.update_tree(self.top_cell)
            progress.Destroy()
                
        dial.Destroy()
        
    def read_seq(self, e):
        adress = self.celltab.get_selected()
        if adress:
            dial = wx.TextEntryDialog(self, 'From which position would you like to translate?', 'Sequence Translation')
            dial.SetValue("0")
            if dial.ShowModal() == wx.ID_OK:
                start = int(dial.GetValue())
                self.top_cell.read_seq(start, energy_per_unit = 0, adress = adress)
                self.celltab.update_tree(self.top_cell)
            dial.Destroy()    
        else:
            dial = wx.MessageDialog(None, 'No cell selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
            
    def split_cell(self, e):
        adress = self.celltab.get_selected()
        if adress != (0,):
            dial = wx.TextEntryDialog(self, 'How much energy should the split cost?', 'Energy Costs')
            dial.SetValue("0")
            if dial.ShowModal() == wx.ID_OK:
                energy = int(dial.GetValue())
                self.top_cell.split_sub_cell(adress, energy)
                self.celltab.update_tree(self.top_cell)
            dial.Destroy()
        else:
            dial = wx.MessageDialog(None, 'Top cell can not be disintegrated', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
    
    def disintegrate_cell(self, e):
        adress = self.celltab.get_selected()
        if adress != (0,):
            dial = wx.TextEntryDialog(self, 'How much energy should the action cost?', 'Energy Costs')
            dial.SetValue("0")
            if dial.ShowModal() == wx.ID_OK:
                energy = int(dial.GetValue())
                self.top_cell.disintegrate_sub_cell(adress, energy)
                self.celltab.update_tree(self.top_cell)
            dial.Destroy()
        else:
            dial = wx.MessageDialog(None, 'Top cell can not be disintegrated', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
            
    def edit_cell_sequence(self, e):
        adress = self.celltab.get_selected()
        if adress:
            sequence = self.top_cell.get_cell_sequence(adress)
            
            dial = wx.TextEntryDialog(self, 'Edit cell sequence:', 
                                      'Edit Cell Sequence', 
                                      style=wx.TE_MULTILINE|wx.OK|wx.CANCEL)
            dial.SetValue(sequence)
            
            if dial.ShowModal() == wx.ID_OK:
                sequence = dial.GetValue()
                self.top_cell.set_cell_sequence(adress, sequence)
                self.celltab.update_tree(self.top_cell)
                
            dial.Destroy()
        else:
            dial = wx.MessageDialog(None, 'No cell selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
            
    def extract_sequence(self, e):
        adress = self.celltab.get_selected()
        if adress:
            sequence = self.top_cell.get_cell_sequence(adress)
            self.add_sequence('', sequence)
        else:
            dial = wx.MessageDialog(None, 'No cell selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
                   
    def save_cell_log(self, e):
        adress = self.celltab.get_selected()
        if adress:
            dlg = wx.FileDialog(self, 'Choose a filename', '', '', '*.csv', wx.OPEN)
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetFilename()
                dirname = dlg.GetDirectory()
                filepath = os.path.join(dirname, filename)
                self.top_cell.save_cell_log_csv(filepath, adress = adress)
                
        else:
            dial = wx.MessageDialog(None, 'No cell selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
    
    def delete_cell(self, e):
        adress = self.celltab.get_selected()
        if adress != False and adress != (0,):
            success = self.top_cell.delete_sub_cell(adress)
            self.celltab.update_tree(self.top_cell)
            if success:
                dial = wx.MessageDialog(None, 'Successfully deleted cell', 'Deleted Cell', wx.OK)
                dial.ShowModal()
            else:
                dial = wx.MessageDialog(None, 'Could not delete cell!', 'Could Not Delete Cell', wx.OK)
                dial.ShowModal()
        elif adress == (0,):
            dial = wx.MessageDialog(None, 'The top cell can not be deleted!\nDelete all cells to achieve the same...', 'Could Not Delete Top Cell', wx.OK)
            dial.ShowModal()
            
        else:
            dial = wx.MessageDialog(None, 'No cell selected...', 'Could Not Perform Action', wx.OK)
            dial.ShowModal()
            
    def flush_all(self, e):
        self.top_cell = cell.Cell(name = 'Top Cell')
        self.celltab.update_tree(self.top_cell)
                            
    def OnQuit(self, e):
        self.Close()
    
    def preferences(self, e):
        adress = self.celltab.get_selected()
        cell = self.top_cell.get_subcell_by_adress(adress)
        
        dial = PreferencesPopup(cell)
        
        if wx.ID_OK == dial.ShowModal():
            sample_mode, sigma_factor, subs_const, interact = dial.get_preferences()
            
            msg = 'Would you like to apply the changes you made:\n\n'
            
            if sample_mode != cell.sample_mode:
                msg += '- Sample Mode = ' + str(sample_mode) + '\n'
            else:
                sample_mode = False
                
            if sigma_factor != cell.normal_factor:
                msg += '- Sigma Factor = ' + str(sigma_factor) + '\n'
            else:
                sigma_factor = False
                
            if subs_const != cell.subs_const:
                msg += '- Constantly added substances:\n'
                for sub, amount in subs_const.items():
                    msg += '\t' + sub + ': ' + str(amount) + ' per round\n'
            else:
                subs_const = False  
                    
            if interact == 'False':
                inter = False
            elif interact == 'True':
                inter = True
                
            if inter != cell.interact_with_subcells:
                msg += '- Cell transport into subcells:  ' + interact + '\n'
            else:
                interact = False
                    
                
            
            msg += '\nto all subcells of the selected cell as well?'
            
            cont = False
            dial = wx.MessageDialog(self, msg, style=wx.YES_NO|wx.CANCEL|wx.ICON_QUESTION)
            res = dial.ShowModal()
            if res == wx.ID_YES:
                cont = True
            dial.Destroy()
            
            #Change parameters in cell
            self.top_cell.change_parameters(adress, sample_mode, sigma_factor, subs_const, cont, interact)
            self.celltab.update_tree(self.top_cell)
                    
    def OnAbout(self, e):
        description = """Sequence Processing 0.1 is a tool to model a pool of cells defined by specific sequences. 
        It entails various methods to run the model and analyze the outcome.
        
        Start by loading an example sequence from a file, creating a cell, 
        translating the sequence in this cell. Then run all cells for some rounds...
        """

        licence = '''Sequence Processing is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        Sequence Processing is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with Sequence Processing.  If not, see <http://www.gnu.org/licenses/>.
        '''

        info = wx.AboutDialogInfo()
        
        info.SetName('Sequence Processing')
        info.SetVersion('0.1')
        info.SetDescription(description)
        info.SetCopyright('(C) 2012 Oliver Dressler')
        info.SetLicence(licence)
        info.AddDeveloper('Oliver Dressler')
        
        wx.AboutBox(info)

def main():
    ex = wx.App()
    SeqGUI(None)
    ex.MainLoop()    

if __name__ == '__main__':
    main()