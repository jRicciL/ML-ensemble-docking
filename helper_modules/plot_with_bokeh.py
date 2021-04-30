import numpy as np
import pandas as pd
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, CDSView, GroupFilter, \
                            Span, CategoricalColorMapper, HoverTool
from bokeh.layouts import row, column
from bokeh.transform import factor_cmap, factor_mark


# Vertical line
vline = Span(location=0, dimension='height', 
             line_color='black', line_width=2, line_alpha=0.5, line_dash='dashed')
# Horizontal line
hline = Span(location=0, dimension='width', 
             line_color='black', line_width=2, line_alpha=0.5, line_dash='dashed')
# HoverTool options
hover= HoverTool(tooltips=[ ('Name', '@name'), ('# Atoms', '@num_atoms'),('Library', '@library'),
          ('Activity', '@Activity')], names = ['actives'])
                            
def create_fig_bokeh(
               desc, 
               source_act, 
               source_inact, 
               col_library_map,
               title='', 
               kind_dr='tsne', 
               legend_location='top_right', 
               legend=False):
    ''' ColumnDataSources source_act  and source_inact must be instantiated'''
    
    f = figure(title=title, plot_width=500, plot_height=420,
          x_axis_label='First Dimension', y_axis_label='Second Dimension',
          tools='pan,box_select,wheel_zoom,reset')
    # Add hovertool 

    f.renderers.extend([vline, hline])
    # Add glyphs
    # Plot inactives
    f_inac = f.circle(x=  f'{desc}_{kind_dr}_x', y=  f'{desc}_{kind_dr}_y', 
               color=col_library_map,
               nonselection_fill_color=col_library_map,
               nonselection_fill_alpha=0.05,
               size=4, alpha=0.3, line_width=0,
               muted_alpha=0.01,
               source=source_inact)

    # Plot actives
    library_names = np.unique(source_act.data['library'])
    df_ = pd.DataFrame(source_act.data)
    for library in library_names:
        data = ColumnDataSource(df_.loc[df_['library'] == library, :])
        f.triangle(x = f'{desc}_{kind_dr}_x', y = f'{desc}_{kind_dr}_y',
               color=col_library_map,
               legend_label=library,
               nonselection_fill_color=col_library_map,
               nonselection_fill_alpha=0.05,
               size=8, line_color='black', line_width=0.5,
               source=data, name=library)
    
    # Styling
    f.title.text_font_size = '1.4em'
    f.axis.axis_label_text_font_size = '1.0em' # font size
    f.axis.axis_label_text_font_style = 'bold'
    f.title.align = 'center'
    f.axis.axis_line_width = 3
    f.axis.major_label_text_font_size = '12pt'
    if legend:
        f.legend.click_policy='hide'
        f.legend.location = legend_location
    else:
        f.legend.visible = False 
    return f