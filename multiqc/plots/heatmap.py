#!/usr/bin/env python
from __future__ import print_function, division, absolute_import

""" MultiQC functions to plot a heatmap """

import json
import logging
import os

from . import get_uid

from multiqc.utils import config

logger = logging.getLogger(__name__)

def plot (data, xcats, ycats=None, pconfig={}):
    """ Plot a 2D heatmap.
    :param data: List of lists, each a representing a row of values.
    :param xcats: Labels for x axis
    :param ycats: Labels for y axis. Defaults to same as x.
    :param pconfig: optional dict with config key:value pairs.
    :return: HTML and JS, ready to be inserted into the page
    """

    if ycats is None:
        ycats = xcats

    # Make a plot
    return highcharts_heatmap(data, xcats, ycats, pconfig)



def highcharts_heatmap (data, xcats, ycats, pconfig={}):
    """
    Build the HTML needed for a HighCharts line graph. Should be
    called by plot_xy_data, which properly formats input data.
    """

    # Reformat the data for highcharts
    pdata = []
    for i, arr in enumerate(data):
        for j, val in enumerate(arr):
            pdata.append([j,i,val])

    # Build the HTML for the page
    if pconfig.get('id') is None:
        pconfig['id'] = 'mqc_hcplot_{}'.format(get_uid())
    html = '<div class="mqc_hcplot_plotgroup">'

    # The 'sort by highlights button'
    html += '''<div class="btn-group hc_switch_group">
        <button type="button" class="mqc_heatmap_sortHighlight btn btn-default btn-sm" data-target="#{id}" disabled="disabled">
            <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
        </button>
        </div>
        '''.format(id=pconfig['id'])

    # The plot div
    html += '''<div class="hc-plot-wrapper">
        <div id="{id}" class="hc-plot not_rendered hc-heatmap"><small>loading..</small></div>
        </div></div>
        '''.format(id=pconfig['id'])

    # Javascript with data dump
    html += '''<script type="text/javascript">
        mqc_plots["{id}"] = {{
            "plot_type": "heatmap",
            "data": {d},
            "xcats": {x},
            "ycats": {y},
            "config": {c}
            }}
        </script>
        '''.format(id=pconfig['id'], d=json.dumps(pdata), x=json.dumps(xcats), y=json.dumps(ycats), c=json.dumps(pconfig));

    return html

