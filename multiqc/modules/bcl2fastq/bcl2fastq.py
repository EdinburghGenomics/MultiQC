from __future__ import division, print_function

from multiqc.modules.base_module import BaseMultiqcModule
import logging
import os
import json
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)

def rat(n, d, nan=None, mul=1.0):
    """ Calculate a ratio while avoiding division by zero errors.
        Strictly speaking we should have nan=float('nan') but for practical
        purposes we'll normally report None (or 0.0?).
    """
    try:
        return ( float(n) * mul ) / float(d)
    except (ZeroDivisionError, TypeError):
        return nan

def pct(n, d, mul=100.0, **kwargs):
    """ Percentage by the same logic.
    """
    return rat(n, d, mul=mul, **kwargs)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcl2fastq', anchor='bcl2fastq',
        href="https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-18.html",
        info="can be used to both demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.")

        self.bcl2fastq_data = dict()

        self.bcl2fastq_bylane = dict()
        self.bcl2fastq_bysample = dict()
        self.source_files = dict()

    def gather(self):
        # Gather data from all json files
        for myfile in self.find_log_files('bcl2fastq'):
            self.parse_file_as_json(myfile)

        # Collect counts by lane and sample (+source_files)
        # This also cleans the sample names as it goes.
        self.split_data_by_lane_and_sample()

        # Filter to strip out ignored sample names
        self.bcl2fastq_bylane = self.ignore_samples(self.bcl2fastq_bylane)
        self.bcl2fastq_bysample = self.ignore_samples(self.bcl2fastq_bysample)

        # Return with Warning if no files are found
        if not ( self.bcl2fastq_bylane or self.bcl2fastq_bysample ):
            log.debug("Could not find any bcl2fastq data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found info for {} lanes, {} samples.".format(
                                 len(self.bcl2fastq_bylane),
                                           len(self.bcl2fastq_bysample) ))

        # Print source files
        for s in self.source_files.keys():
            self.add_data_source(s_name=s, source=",".join(list(set(self.source_files[s]))), module='bcl2fastq', section='bcl2fastq-bysample')

        # Add sample counts to general stats table
        self.add_general_stats()
        self.write_data_file({str(k): self.bcl2fastq_bylane[k] for k in self.bcl2fastq_bylane.keys()}, 'multiqc_bcl2fastq_bylane')
        self.write_data_file(self.bcl2fastq_bysample, 'multiqc_bcl2fastq_bysample')

        # Add section for summary stats per flow cell
        self.add_section (
            name = 'Lane Statistics',
            anchor = 'bcl2fastq-lanestats',
            description = 'Statistics about each lane for each flowcell',
            plot = self.lane_stats_table()
        )

        # Add section for counts by lane
        cats = OrderedDict()
        cats["perfect"] = {'name': 'Perfect Index Reads'}
        cats["imperfect"] = {'name': 'Mismatched Index Reads'}
        cats["undetermined"] = {'name': 'Undetermined Reads'}
        self.add_section (
            name = 'Clusters by lane',
            anchor = 'bcl2fastq-bylane',
            description = 'Number of reads per lane (with number of perfect index reads)',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.
                To avoid conflicts the run ID is prepended.""",
            plot = bargraph.plot(self.get_bar_data_from_counts(self.bcl2fastq_bylane), cats)
        )

        # Add section for counts by sample
        self.add_section (
            name = 'Clusters by sample',
            anchor = 'bcl2fastq-bysample',
            description = 'Number of reads per sample (with number of perfect index reads)',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combinned. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.
                To avoid conflicts the runId is prepended.""",
            plot = bargraph.plot(self.get_bar_data_from_counts(self.bcl2fastq_bysample), cats)
        )

    def parse_file_as_json(self, myfile):
        try:
            content = json.loads(myfile["f"])
        except ValueError:
            log.warn('Could not parse file as json: {}'.format(myfile["fn"]))
            return
        runId = content["RunId"]
        if not runId in self.bcl2fastq_data:
            self.bcl2fastq_data[runId] = dict()
        run_data = self.bcl2fastq_data[runId]
        for conversionResult in content["ConversionResults"]:
            lane = 'L{}'.format(conversionResult["LaneNumber"])
            if lane in run_data:
                log.debug("Duplicate runId/lane combination found! Overwriting: {}".format(self.prepend_runid(runId, lane)))
            run_data[lane] = {
                "total": 0,
                "total_yield": 0,
                "perfectIndex": 0,
                "samples": dict(),
                "yieldQ30": 0,
                "qscore_sum": 0
            }
            for demuxResult in conversionResult["DemuxResults"]:
                sample = demuxResult["SampleName"]
                if sample in run_data[lane]["samples"]:
                    log.debug("Duplicate runId/lane/sample combination found! Overwriting: {}, {}".format(self.prepend_runid(runId, lane),sample))
                run_data[lane]["samples"][sample] = {
                    "total": 0,
                    "total_yield": 0,
                    "perfectIndex": 0,
                    "filename": os.path.join(myfile['root'],myfile["fn"]),
                    "yieldQ30": 0,
                    "qscore_sum": 0
                }
                run_data[lane]["total"] += demuxResult["NumberReads"]
                run_data[lane]["total_yield"] += demuxResult["Yield"]
                run_data[lane]["samples"][sample]["total"] += demuxResult["NumberReads"]
                run_data[lane]["samples"][sample]["total_yield"] += demuxResult["Yield"]
                #If there was no index, there will be no metrics!
                indexMetric = None
                for indexMetric in demuxResult.get("IndexMetrics",[]):
                    run_data[lane]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                    run_data[lane]["samples"][sample]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                    run_data[lane]["samples"][sample]["IndexSequence"] = indexMetric["IndexSequence"]
                if indexMetric is None:
                    run_data[lane]["perfectIndex"] = None
                    run_data[lane]["samples"][sample]["perfectIndex"] = None
                for readMetric in demuxResult["ReadMetrics"]:
                    run_data[lane]["yieldQ30"] += readMetric["YieldQ30"]
                    run_data[lane]["qscore_sum"] += readMetric["QualityScoreSum"]
                    run_data[lane]["samples"][sample]["yieldQ30"] += readMetric["YieldQ30"]
                    run_data[lane]["samples"][sample]["qscore_sum"] += readMetric["QualityScoreSum"]
            undeterminedYieldQ30 = 0
            undeterminedQscoreSum = 0
            #If there was no index, there will be no undetermined metrics.
            conversionResult.setdefault("Undetermined", dict( ReadMetrics = [],
                                                              NumberReads = 0,
                                                              Yield = 0.0 ))

            for readMetric in conversionResult["Undetermined"]["ReadMetrics"]:
                undeterminedYieldQ30 += readMetric["YieldQ30"]
                undeterminedQscoreSum += readMetric["QualityScoreSum"]

            run_data[lane]["samples"]["undetermined"] = {
                "total": conversionResult["Undetermined"]["NumberReads"],
                "total_yield": conversionResult["Undetermined"]["Yield"],
                "perfectIndex": 0,
                "yieldQ30": undeterminedYieldQ30,
                "qscore_sum": undeterminedQscoreSum
            }

        # Calculate Percents and averages. Beware of division by zero!
        for lane in run_data:
            run_data[lane]["percent_Q30"] = pct(run_data[lane]["yieldQ30"], run_data[lane]["total_yield"])
            run_data[lane]["percent_perfectIndex"] = pct(run_data[lane]["perfectIndex"], run_data[lane]["total"])
            run_data[lane]["mean_qscore"] = rat(run_data[lane]["qscore_sum"], run_data[lane]["total_yield"])
            for sample, d in run_data[lane]["samples"].items():
                run_data[lane]["samples"][sample]["percent_Q30"] = pct(d["yieldQ30"], d["total_yield"])
                run_data[lane]["samples"][sample]["percent_perfectIndex"] = pct(d["perfectIndex"], d["total"])
                run_data[lane]["samples"][sample]["mean_qscore"] = rat(d["qscore_sum"], d["total_yield"])

    def split_data_by_lane_and_sample(self):
        for runId in self.bcl2fastq_data.keys():
            for lane in self.bcl2fastq_data[runId].keys():
                uniqLaneName = self.prepend_runid(runId, lane)
                self.bcl2fastq_bylane[uniqLaneName] = {
                    "total": self.bcl2fastq_data[runId][lane]["total"],
                    "total_yield": self.bcl2fastq_data[runId][lane]["total_yield"],
                    "perfectIndex": self.bcl2fastq_data[runId][lane]["perfectIndex"],
                    "undetermined": self.bcl2fastq_data[runId][lane]["samples"]["undetermined"]["total"],
                    "yieldQ30": self.bcl2fastq_data[runId][lane]["yieldQ30"],
                    "qscore_sum": self.bcl2fastq_data[runId][lane]["qscore_sum"],
                    "percent_Q30": self.bcl2fastq_data[runId][lane]["percent_Q30"],
                    "percent_perfectIndex": self.bcl2fastq_data[runId][lane]["percent_perfectIndex"],
                    "mean_qscore": self.bcl2fastq_data[runId][lane]["mean_qscore"]
                }

                for sample, sinfo in self.bcl2fastq_data[runId][lane]["samples"].items():
                    # Note - at this point the sample names are uncleaned.
                    # so clean it up
                    csample = self.clean_s_name(sample, '.')

                    totals = self.bcl2fastq_bysample.setdefault(csample, {
                            "total": 0,
                            "total_yield": 0,
                            "perfectIndex": None,
                            "yieldQ30": 0,
                            "qscore_sum": 0,
                            "IndexSequence": set(),
                            "SampleNames" : list(),
                        })

                    totals["total"] += sinfo["total"]
                    totals["total_yield"] += sinfo["total_yield"]
                    if sinfo["perfectIndex"] is not None:
                        totals["perfectIndex"] = (totals["perfectIndex"] or 0) + sinfo["perfectIndex"]
                    totals["yieldQ30"] += sinfo["yieldQ30"]
                    totals["qscore_sum"] += sinfo["qscore_sum"]
                    if sinfo.get("IndexSequence"):
                        totals["IndexSequence"].add(sinfo["IndexSequence"])
                    # Add the cleaned sample name so I can infer the pool. I may need to explicitly
                    # capture SapleID instead if we change the format of the sample sheets at all.
                    totals["SampleNames"].append(sample)

                    totals["percent_Q30"]          = pct(totals["yieldQ30"], totals["total_yield"])
                    totals["percent_perfectIndex"] = pct(totals["perfectIndex"], totals["total"])
                    totals["mean_qscore"]          = rat(totals["qscore_sum"], totals["total_yield"])
                    if sinfo.get("filename"):
                        self.source_files.setdefault(csample,[]).append( sinfo["filename"] )

    def add_general_stats(self):
        data = {
            key: {
                "yieldQ30": val["yieldQ30"],
                "total": val["total"],
                "perfectPercent": pct( val["perfectIndex"], val["total"] )
            } for key, val in self.bcl2fastq_bysample.items()
        }
        headers = OrderedDict()
        headers['total'] = {
            'title': 'Total Fragments',
            'description': 'Total number of read pairs (or single-end reads)' + \
                           'for this sample as determined by bcl2fastq demultiplexing ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'Blues',
            'shared_key': 'read_count',
            'format': '{:,}'
        }

        # I'm hoping I can excise this but I've been asked to try including it...
        # We could show that as a percentage?  I have to calculate these explicitly.
        # Also this may not make sense outside of our assumption that we are processing one lane at a time.
        grand_total_reads = sum( s["total"] for s in self.bcl2fastq_bysample.values() )

        # This should always be the same, right???
        grand_total_2 = sum( lane["total"] for run in self.bcl2fastq_data.values() for lane in run.values() )
        assert grand_total_2 == grand_total_reads

        for key, val in self.bcl2fastq_bysample.items():
            data[key]['total_as_pct'] = pct( val["total"], grand_total_reads )

        headers['total_as_pct'] = {
            'title': '% Fragments',
            'description': 'Total number of fragments expressed as a percentage of the total processed.',
            'max': 100,
            'min': 0,
            'scale': 'Blues',
            'format': '{0:.2f}'
        }
        # ...End of dubious bit

        headers['yieldQ30'] = {
            'title': '{} Yield &ge; Q30'.format(config.base_count_prefix),
            'description': 'Number of bases with a Phred score of 30 or higher ({})'.format(config.base_count_desc),
            'min': 0,
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count'
        }
        headers['perfectPercent'] = {
            'title': '% Perfect Index',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%',
            'format': '{0:.1f}'
        }

        # Now add the index sequences to the data. In the general case, it's possible that the
        # sample appeared on one lane with multiple indexes or on multiple lanes with different
        # indexes so we just join the set with commas.
        if getattr(config, 'bcl2fastq_config', {}).get('add_index_sequences'):
            for k, v in data.items():
                v['indexSeq'] = ','.join(sorted(self.bcl2fastq_bysample[k].get("IndexSequence", [])))

        if any(v.get('indexSeq') for v in data.values()):
            headers['indexSeq'] = {
                'title': 'Index Sequence',
                'description': 'Index sequence as set in the Sample Sheet (ie. as read by the sequencer)',
                'format': '{s}',
                'placement': 5.0,
                'textcell': True
            }

        # And likewise for the pool name (yes, this is very Edinburgh Genomics centric but I need it
        # to work so I'm stuffing it in here!)
        # Again, it's theoretically possible to see the same sample in two pools so use a set.
        if getattr(config, 'bcl2fastq_config', {}).get('add_pool_names'):
            for k, v in data.items():
                pools_for_sample = [ sn.split("__")[0] if "__" in sn else ""
                                     for sn in self.bcl2fastq_bysample[k].get("SampleNames", []) ]
                pools_set = set( pn for pn in pools_for_sample
                                 if pn.lower not in ['', 'none', 'nopool'] )
                v['pool'] = ','.join(sorted(pools_set))

        if any(v.get('pool') for v in data.values()):
            headers['pool'] = {
                'title': 'Pool Name',
                'description': 'Pool in which this library was sequenced.',
                'format': '{s}',
                'placement': 4.0,
                'textcell': True
            }

        self.general_stats_addcols(data, headers)

    def lane_stats_table(self):
        """ Return a table with overview stats for each bcl2fastq lane for a single flow cell """
        headers = OrderedDict()
        headers['total_yield'] = {
            'title': '{} Total Yield'.format(config.base_count_prefix),
            'description': 'Number of bases ({})'.format(config.base_count_desc),
            'min': 0,
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count'
        }
        headers['total'] = {
            'title': '{} Total Clusters'.format(config.read_count_prefix),
            'description': 'Total number of clusters for this lane ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        headers['percent_Q30'] = {
            'title': '% bases &ge; Q30',
            'description': 'Percentage of bases with greater than or equal to Q30 quality score',
            'suffix': '%',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn'
        }
        headers['mean_qscore'] = {
            'title': 'Mean Quality',
            'description': 'Average phred qualty score',
            'min': 0,
            'scale': 'Spectral'
        }
        headers['percent_perfectIndex'] = {
            'title': '% Perfect Index',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        table_config = {
            'namespace': 'bcl2fastq',
            'id': 'bcl2fastq-lane-stats-table',
            'table_title': 'bcl2fastq Lane Statistics',
            'col1_header': 'Run ID - Lane',
            'no_beeswarm': True
        }
        return table.plot(self.bcl2fastq_bylane, headers, table_config)

    def prepend_runid(self, runId, rest):
        return str(runId)+" - "+str(rest)

    def get_bar_data_from_counts(self, counts):
        bar_data = {}
        for key, value in counts.items():
            if value["perfectIndex"] is not None:
                bar_data[key] = {
                    "perfect": value["perfectIndex"],
                    "imperfect": value["total"] - value["perfectIndex"],
                }
            else:
                bar_data[key] = {
                    "noindex": value["total"]
                }
            if "undetermined" in value:
                bar_data[key]["undetermined"] = value["undetermined"]
        return bar_data
