// Listen for clicks on table cells
$(document).on('click', '#alignment_table table.dataTable tbody td', function() {
    var table = $('#alignment_table table.dataTable').DataTable();
    var cell = table.cell(this);
    var colIndex = cell.index().column;
    var rowIndex = cell.index().row;
    
    // Remove previous column highlights
    $('#alignment_table table.dataTable tbody td').removeClass('selected-column');

    // Highlight all cells in the selected column
    $('#alignment_table table.dataTable tbody tr').each(function() {
        $(this).find('td').eq(colIndex).addClass('selected-column');
    });

    // Remove previous row highlights
    $('#alignment_table table.dataTable tbody tr td').removeClass('selected-row');

    // Highlight the clicked row
    var rowElement = $(this).closest('tr');
    rowElement.find('td').addClass('selected-row');

    // Remove previous cell border
    $('#alignment_table table.dataTable tbody td').removeClass('selected-cell');

    // Add border to the selected cell
    $(this).addClass('selected-cell');

    // Send column and row indices to Shiny
    Shiny.setInputValue('col_selected', colIndex, {priority: 'event'});
    Shiny.setInputValue('row_selected', rowIndex, {priority: 'event'});
});

// Handler for scrolling table to specific column
Shiny.addCustomMessageHandler('scrollTableToColumn', function(columnIndex) {
    console.log('scrollTableToColumn message received with columnIndex:', columnIndex);

    // Find the DataTable's scrollable container and table element
    var scrollContainer = document.querySelector('#alignment_table .dataTables_scrollBody');
    var table = document.querySelector('#alignment_table table.dataTable');

    if (scrollContainer && table) {
        console.log('✅ Found table');

        // Use DataTables API to get rows
        var dt = $('#alignment_table table.dataTable').DataTable();
        var rowNodes = dt.rows().nodes();
        
        // Get columns from the first row instead of the header
        var columns = rowNodes.length > 0 ? rowNodes[0].querySelectorAll('td') : [];
        
        console.log('✅ Column count based on first row:', columns.length);

        // Calculate scroll position and perform scroll if column index is valid
        if (columnIndex < columns.length) {
            var offset = 0;
            
            // Sum widths of preceding columns to get scroll position
            for (var i = 0; i < columnIndex; i++) {
                offset += columns[i].offsetWidth;
            }
            
            console.log('✅ Calculated offset for column:', offset);
            scrollContainer.scrollLeft = offset;
        } else {
            console.log('❌ Invalid columnIndex:', columnIndex);
        }
    } else {
        console.log('❌ Scrollable container or table not found!');
    }
});
