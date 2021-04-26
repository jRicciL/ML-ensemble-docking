# Libraries
import numpy as np 
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt
import pylab
import warnings

# Metric imported from other libraries
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc

class PlotMetric:
    def __init__(self, y_true, df_predictions, decreasing = True, color_palette = 'Dark2', figsize = (7,7)):
        if type(y_true) is not np.ndarray:
            try: 
               self.y_true = np.array(y_true)
            except ValueError:
                print('y_true should be a numpy array with values 1 = active and 0 = inactive')
        else:
            self.y_true = y_true
        if not np.array_equal(y_true, y_true.astype(bool)):
            assert 'y_true array must be binary'
        if type(df_predictions) is not pd.DataFrame or len(df_predictions) < 1:
             raise ValueError('Error_ temp text')

        
        self.N = len(y_true)
        self.n = len(y_true[y_true == 1])
        self.R_a = self.n/self.N
        self.df_predictions = df_predictions.copy()

        # We make sure that preds are numpy arrays
#         for key, y_pred in self.df_predictions.iteritems():
#             if decreasing:
#                 self.df_predictions[key] = -1 * np.array(y_pred)
#             else:
#                 self.df_predictions[key] = np.array(y_pred)
        if decreasing:
            self.df_predictions =+ self.df_predictions*(-1)
        
        self.color_palette = color_palette
        self.available_metrics = {'roc_auc': self._get_roc_auc,
                                  'p_roc': self._get_pRoc_auc,
                                  'auac': self._get_ac_auc,
                                  'pr_auc': self._get_pr_auc,
                                  'ref_auc': self._get_ref_auc,
                                  'nef_auc': self._get_nef_auc,
                                  'ef': self._get_ef_value,
                                  'rie': self._rie_score,
                                  'bedroc': self._bedroc_score}
        pylab.rcParams['figure.figsize'] = figsize
        sns.set( context = 'talk', style = 'white', palette = color_palette)

    def get_metrics_report(self, omit = None):
        pass
        #if type(omit) == str and omit in self.available_metrics
        # TODO
    
    # ROC
    def _get_roc(self, y_pred):
        N = self.N
        n = self.n
        order = np.argsort(- y_pred)
        y_true_ord = self.y_true[order]
        # Claculo de  TPR
        tpr = np.zeros(N)
        fpr = np.zeros(N)
        n_actives = 0
        n_inactives = 0
        for i, active in enumerate(y_true_ord):
            if active:
                n_actives += 1
            else:
                n_inactives += 1
            tpr[i] = n_actives
            fpr[i] = n_inactives
        # Normalizamos los valores
        if n_actives > 0:
            tpr = tpr/n_actives
        if n_inactives > 0:
            fpr = fpr/n_inactives
        # The we insert 0 at the begining of the list
        fpr = np.insert(fpr, 0, 0)
        tpr = np.insert(tpr, 0, 0)
        #fpr, tpr, _ = roc_curve(self.y_true, y_pred)
        return fpr, tpr
    
    # ROC-AUC
    def _get_roc_auc(self, y_pred):
        return(roc_auc_score(y_true = self.y_true, y_score = y_pred))

    # pROC
    def _get_pRoc(self, y_pred):
        fpr, tpr = self._get_roc(y_pred = y_pred)
        # Clark correction suggest to change fpr(i) by (1/N) if fpr(i) = 0
        p_fpr = [(- np.log10(1/i)) if i != 0 else - np.log10(self.N) for i in fpr]
        return p_fpr, tpr

    def _get_pRoc_auc(self, y_pred, normalized=True):
        p_fpr, tpr = self._get_pRoc(y_pred = y_pred)
        # I suggest to normalize the value by using the number of molecules
        pRoc = auc(p_fpr, tpr)
        
        if normalized:
            pRoc = pRoc / (- np.log10(1 / self.N))
        return pRoc

    # Accumulation Curve
    def _get_ac(self, y_pred, normalized):
        N = self.N
        n = self.n
        order = np.argsort(- y_pred)
        y_true_ord = self.y_true[order]
        ranking_pos = np.linspace(0, N, N + 1)
        n_counter = 0
        f_k = np.zeros(N + 1)
        for i, k in enumerate(y_true_ord):
            if k: # If is active
                n_counter += 1
            f_k[i + 1] = n_counter
        if normalized:
            ranking_pos = ranking_pos / N
            f_k = f_k / n
        return ranking_pos, f_k

    def _get_ac_auc(self, y_pred, normalized = True):
        k, f_k = self._get_ac(y_pred = y_pred, normalized = normalized)
        auac = auc(k, f_k)
        return auac

    # Precision and Recall
    def _get_pr(self, y_pred):
        precision, recall, thresholds = precision_recall_curve(y_true = self.y_true, probas_pred = y_pred)
        return precision, recall, thresholds
    # PR-AUC
    def _get_pr_auc(self, y_pred):
        precision, recall, thresholds = self._get_pr(y_pred = y_pred)
        pr_auc = auc(recall, precision)
        return pr_auc

    # Enrichment Factor
    def _get_ef(self, y_pred, fractions, method='normalized'):
        N = self.N
        n = self.n
        order = np.argsort(- y_pred)
        y_true_ord = self.y_true[order]
        
        efs = []
        if type(fractions) is np.ndarray:
            fractions = fractions.tolist()
        # Fractions mus be ordered
        fractions.sort()
        # Get the N_s values
        N_s_floor = [np.floor(N * f) for f in fractions]
        # If the fraction 1.0 is not in fractions we need to add it
        if fractions[-1] != 1:
            N_s_floor.append(N)
        # if 0 is given in the fractions:
        if N_s_floor[0] == 0:
            efs.append(0)
            N_s_floor.pop(0)
        # start the counting of actives at 0
        n_s = 0
        for n_mol in range(N):
            if n_mol > N_s_floor[0] and n_mol > 0:
                N_s = n_mol
                if method == 'relative':
                    ef_i = (100 * n_s) / min(N_s, n)
                elif method == 'absolute':
                    ef_i = (N * n_s) / (n * N_s)
                else:
                    ef_i = n_s / min(N_s, n)
                efs.append(ef_i)
                N_s_floor.pop(0)
            active = y_true_ord[n_mol]
            if active:
                n_s += 1
        if fractions[-1] == 1 and N_s_floor[0] == N:
            if method == 'relative': 
                ef_i = 100
            else: 
                ef_i = 1
            efs.append(ef_i)
        return efs
    
    def get_efs(self, method, fractions = [0.005, 0.01, 0.02, 0.05], rounded = 2):
        ef_results = {}
        for key, y_pred in self.df_predictions.iteritems():
            ef_results[key] = self._get_ef(y_pred, method = method, fractions = fractions)
        names = {'relative': 'REF', 'absolute': 'EF', 'normalized': 'NEF'}
        row_names_ef = [F'{names[method]} at {i*100}%' for i in fractions]
        df_efs = pd.DataFrame(ef_results, index = row_names_ef)
        df_efs["#ligs at X%"] = [np.floor(i* self.N) for i in fractions]
        df_efs = df_efs.round(rounded)
        return df_efs.T

    def _get_ef_value(self, y_pred, fraction, method='absolute'):
        ef_value = self._get_ef(y_pred = y_pred, method = method, fractions = [fraction])
        return ef_value

    def _get_ref_auc(self, y_pred, method = 'relative', max_chi = 1):
        fractions = np.linspace(0.0, max_chi, len(self.y_true) - 2 )
        efs = self._get_ef(y_pred = y_pred, method = method, fractions = fractions)
        efs_auc = auc(fractions, efs)
        return efs_auc

    def _get_nef_auc(self, y_pred, method = 'normalized', max_chi = 1):
        fractions = np.linspace(0.0, max_chi, len(self.y_true) - 2 )
        efs = self._get_ef(y_pred = y_pred, method = method, fractions = fractions)
        efs_auc = auc(fractions, efs)
        return efs_auc

    def _rie_score(self, y_pred, alpha):
        N = self.N
        n = self.n
        y_true = self.y_true
        
        warnings.simplefilter('always', UserWarning)
        if alpha*self.R_a >= 1:
            max_alpha = 1/self.R_a
            warnings.warn('Parsed alpha value ({:0.2f}) times R_a is greater than 1.\nAn alpha value below {:0.2f} is recommended.'.format(alpha, max_alpha),
            ResourceWarning)
        order = np.argsort(- y_pred)
        m_rank = (y_true[order] == 1).nonzero()[0] + 1
        S_numerator = np.sum(np.exp(-alpha * m_rank / N))
        _denominator =  (n / N) * ( (1 - np.exp(-alpha)) / (np.exp(alpha/N) - 1) )
        
        RIE = S_numerator / _denominator
        return RIE

    def _rie_max_min_scores(self, alpha):
        R_a = self.R_a
        RIE_min = (1 - np.exp(alpha * R_a)) / (R_a*(1 - np.exp(alpha)))
        RIE_max = (1 - np.exp(- alpha * R_a)) / (R_a*(1 - np.exp(- alpha)))
        return (RIE_min, RIE_max)

    def _bedroc_score(self, y_pred, alpha):
        N = self.N
        n = self.n
        y_true = self.y_true
        R_a = self.R_a
        
        RIE = self._rie_score(y_pred, alpha)
        RIE_min, RIE_max = self._rie_max_min_scores(alpha)
        
        BEDROC = (RIE - RIE_min) / (RIE_max - RIE_min)
        return BEDROC

    def get_ries(self, alphas = [5, 10, 20], rounded = 3):
        rie_results = {}
        for key, y_pred in self.df_predictions.iteritems():
            rie_results[key] = [self._rie_score(y_pred = y_pred, alpha = a) for a in alphas]
        row_names_rie = [F'alpha = {a}' for a in alphas]
        df_rie = pd.DataFrame(rie_results, index = row_names_rie)
        df_rie = df_rie.round(rounded).T
        return df_rie

    def get_bedrocs(self, alphas = [5, 10, 20], rounded = 3):
        bedrocs_results = {}
        for key, y_pred in self.df_predictions.iteritems():
            bedrocs_results[key] = [self._bedroc_score(y_pred = y_pred, alpha = a) for a in alphas]
        row_names_bedroc = [F'alpha = {a}' for a in alphas]
        df_bedrocs = pd.DataFrame(bedrocs_results, index = row_names_bedroc)
        df_bedrocs = df_bedrocs.round(rounded).T
        return df_bedrocs

    ## PLOTTING FUNCTIONS
    def _add_plot_ef(self, y_pred, label, max_chi, max_num_of_ligands, method = 'normalized',  **kwargs):
        fractions = np.linspace(0.0, 1, len(self.y_true) - 2 )
        names = {'relative': 'REF', 'absolute': 'EF', 'normalized': 'NEF'}
        efs = self._get_ef(y_pred = y_pred, method = method, fractions = fractions)
        # Calculate max # of ranked values to show
        if max_num_of_ligands is not None and max_chi == 1:
            n_rankers = max_num_of_ligands
            fractions = range(n_rankers)
            efs = efs[:n_rankers]
        else:
            n_rankers = int(np.ceil(max_chi*len(fractions)))
            fractions = fractions[:n_rankers]
            efs =  efs[:n_rankers]
        efs_auc = auc(fractions, efs)
        if method == 'absolute':
            label = label + F' {names[method]}'
        else:
            label = label + F' {names[method]}-AUC' + ' = %0.2f' % efs_auc
        plt.plot(fractions, efs, label = label, **kwargs)

    def plot_ef_auc(self, title, method = 'normalized', 
                    max_chi = 1, max_num_of_ligands = None,
                    keys_to_omit = [], key_to_plot = None, key_to_fade = None,
                     fontsize='x-small', showplot = True, show_by_itself=True, **kwargs):

        methods = ('relative', 'absolute', 'normalized')
        method = method.lower()
        if method not in methods:
            raise AttributeError(F'method value, {method} is not available.\nAvailable methods are:\n{methods}')
        sns.color_palette(self.color_palette)

        if method == 'absolute' and max_num_of_ligands is not None:
            raise AttributeError(F'arguments method="absolute" and "max_num_of_ligands" can not be applied together.')
        for key, y_pred in self.df_predictions.iteritems():
            if key in keys_to_omit:
                continue
            if key_to_fade == key:
                self._add_plot_ef(y_pred, method = method, label = str(key), 
                        max_chi = max_chi, max_num_of_ligands = max_num_of_ligands,
                        linestyle = '--', 
                        linewidth = 1.5, **kwargs)
                continue
            if type(key_to_plot) is str and key_to_plot in self.df_predictions.keys():
                key = key_to_plot
                y_pred = self.df_predictions[key]
                self._add_plot_ef(y_pred, method = method, label = str(key), 
                                    max_chi = max_chi, max_num_of_ligands = max_num_of_ligands, **kwargs)
                break
            self._add_plot_ef(y_pred, method = method, label = str(key), 
                                    max_chi = max_chi, max_num_of_ligands = max_num_of_ligands, **kwargs)
        if showplot:
            plt.legend(fontsize=fontsize)
            if method == 'absolute':
                plt.plot([self.R_a, self.R_a], [0 , 1/self.R_a], 'k--', c = 'grey')
                #plt.plot([0, 1], [1, 1], 'k--', c = 'grey')
            else: plt.ylim(0, 1.1)
            if max_num_of_ligands is not None and max_chi == 1:
                plt.xlabel("# ligands at ranking top")
            else:
                plt.xlabel("Ranking Fraction")
            plt.ylabel(F"Enrichment Factor ({method[:3].capitalize()}.)")
            
            plt.grid(linestyle='--', linewidth='0.8')
            plt.title(title)
            if show_by_itself:
                plt.show()

    def _add_plot_pr(self, y_pred, label, **kwargs):
        precision, recall, thresholds = self._get_pr(y_pred)
        auc_pr = self._get_pr_auc(y_pred)
        plt.plot(recall, precision, label = label + ' AUC-PR = %0.3f' % auc_pr, **kwargs)

    def plot_pr_auc(self, title, keys_to_omit = [], key_to_plot = None, key_to_fade = None,
                     fontsize='x-small', showplot = True, show_by_itself=True, **kwargs):
        sns.color_palette(self.color_palette)
        
        for key, y_pred in self.df_predictions.iteritems():
            if key in keys_to_omit:
                continue
            if key_to_fade == key:
                self._add_plot_pr(y_pred, label = str(key), linestyle = '--', 
                        linewidth = 1.5, **kwargs)
                continue
            if type(key_to_plot) is str and key_to_plot in self.df_predictions.keys():
                key = key_to_plot
                y_pred = self.df_predictions[key]
                self._add_plot_pr(y_pred, label = str(key), **kwargs)
                break
            self._add_plot_pr(y_pred, label = str(key), **kwargs)
        if showplot:
            plt.legend(fontsize=fontsize)
            no_skill = self.R_a
            plt.plot([0, 1], [no_skill, no_skill], 'k--')
            plt.xlabel("Recall")
            plt.ylabel("Precision")
            plt.ylim(0, 1.1)
            plt.grid(linestyle='--', linewidth='0.8')
            plt.title(title)
            if show_by_itself:
                plt.show()

    # pROC
    def _add_plot_pRoc(self, y_pred, label, **kwargs):
        fpr, tpr = self._get_pRoc(y_pred)
        auc = self._get_pRoc_auc(y_pred)
        plt.plot(fpr, tpr, label = label + ' AUC-pROC = %0.3f' % auc, **kwargs)
        
    def plot_pRoc_auc(self, title, keys_to_omit = [], key_to_plot = None, key_to_fade = None,
                     fontsize='small', showplot = True, show_by_itself=True, **kwargs):
        sns.color_palette(self.color_palette)
        
        for key, y_pred in self.df_predictions.iteritems():
            if key in keys_to_omit:
                continue
            if key_to_fade == key:
                self._add_plot_pRoc(y_pred, label = str(key), linestyle = '--', 
                        linewidth = 1.5, **kwargs)
                continue
            if type(key_to_plot) is str and key_to_plot in self.df_predictions.keys():
                key = key_to_plot
                y_pred = self.df_predictions[key]
                self._add_plot_pRoc(y_pred, label = str(key), **kwargs)
                break
            self._add_plot_pRoc(y_pred, label = str(key), **kwargs)
        if showplot:
            plt.legend(fontsize=fontsize)
            #plt.plot([0, 1], [0, 1], 'k--', c = 'gray')
            plt.xlabel("pFPR (log10(1 - specificity))")
            plt.ylabel("TPR (sensitivity)")
            plt.grid(linestyle='--', linewidth='0.8')
            plt.title(title)
            if show_by_itself:
                plt.show()

    def _add_plot_bedroc(self, label, alphas, bedrocs, **kwargs):
        plt.plot(alphas, bedrocs, label = label,  **kwargs)

    def plot_bedroc(self, title, alphas, keys_to_omit = [], key_to_plot = None, key_to_fade = None,
                     fontsize='small', showplot = True, show_by_itself=True, **kwargs):
        sns.color_palette(self.color_palette)
        df_bedroc = self.get_bedrocs(alphas = alphas).T
        for key in self.df_predictions.keys():
            if key in keys_to_omit:
                continue
            if key_to_fade == key:
                self._add_plot_bedroc(label = str(key), alphas = alphas, linestyle = '--', marker = 'o',
                bedrocs = df_bedroc[key].values, **kwargs)
                continue
            if type(key_to_plot) is str and key_to_plot in self.df_predictions.keys():
                key = key_to_plot
                self._add_plot_bedroc(label = str(key), alphas = alphas, bedrocs = df_bedroc[key].values, marker = 'v', **kwargs)
                break
            self._add_plot_bedroc(label = str(key), alphas = alphas, bedrocs = df_bedroc[key].values, marker = 'v', **kwargs)
        if showplot:
            plt.legend(fontsize=fontsize)
            plt.xlabel("Alpha values")
            plt.ylabel("BEDROC score")
            plt.grid(linestyle='--', linewidth='0.8')
            plt.title(title)
            plt.ylim(-0.1, 1.1)
            if show_by_itself:
                plt.show()

    # Plotting functions
    def _add_plot_roc(self, y_pred, label, **kwargs):
        fpr, tpr = self._get_roc(y_pred)
        auc = self._get_roc_auc(y_pred)
        plt.plot(fpr, tpr, label = label + ' AUC-ROC = %0.3f' % auc, **kwargs)
        
    def plot_roc_auc(self, title, keys_to_omit = [], key_to_plot = None, key_to_fade = None,
                     fontsize='small', showplot = True, show_by_itself=True, **kwargs):
        sns.color_palette(self.color_palette)
        
        for key, y_pred in self.df_predictions.iteritems():
            if key in keys_to_omit:
                continue
            if key_to_fade == key:
                self._add_plot_roc(y_pred, label = str(key), linestyle = '--', 
                        linewidth = 1.5, **kwargs)
                continue
            if type(key_to_plot) is str and key_to_plot in self.df_predictions.keys():
                key = key_to_plot
                y_pred = self.df_predictions[key]
                self._add_plot_roc(y_pred, label = str(key), **kwargs)
                break
            self._add_plot_roc(y_pred, label = str(key), **kwargs)
        if showplot:
            plt.legend(fontsize=fontsize)
            plt.plot([0, 1], [0, 1], 'k--', c = 'gray')
            plt.xlabel("FPR (1 - specificity)")
            plt.ylabel("TPR (sensitivity)")
            plt.grid(linestyle='--', linewidth='0.8')
            plt.title(title)
            if show_by_itself:
                plt.show()

    def _add_plot_ac(self, y_pred, label, normalized = True, **kwargs):
        k, f_k = self._get_ac(y_pred = y_pred, normalized = normalized)
        auc_ac = self._get_ac_auc(y_pred = y_pred, normalized = normalized)
        plt.plot(k, f_k , label = label + ' AUAC = %0.3f' % auc_ac, **kwargs)

    def plot_auac(self, title, keys_to_omit = [], key_to_plot = None, key_to_fade = None, normalized = True,
                     fontsize='small', showplot = True, show_by_itself=True, **kwargs):
        sns.color_palette(self.color_palette)
        
        for key, y_pred in self.df_predictions.iteritems():
            if key in keys_to_omit:
                continue
            if key_to_fade == key:
                self._add_plot_ac(y_pred, label = str(key), normalized = normalized,
                linestyle = '--', linewidth = 1.5, **kwargs)
                continue
            if type(key_to_plot) is str and key_to_plot in self.df_predictions.keys():
                key = key_to_plot
                y_pred = self.df_predictions[key]
                self._add_plot_ac(y_pred, label = str(key), normalized = normalized, **kwargs)
                break
            self._add_plot_ac(y_pred, label = str(key), normalized = normalized, **kwargs)
        if showplot:
            plt.legend(fontsize=fontsize)
            if normalized:
                x_label = 'Normalized ranking positions (x = k/N)'
                y_label = 'Normalized num. of actives (Fa(x))'
            else:
                x_label = 'Ranking positions (k)'
                y_label = 'Num. of actives (Fa(k))'
            plt.plot([0, 0], [1, 1], 'k--', c = 'grey')
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.grid(linestyle='--', linewidth='0.8')
            plt.title(title)
            if show_by_itself:
                plt.show()
        
    # Plotting distributions
    def plot_actives_distribution(self, colors = {1: '#e74c3c', 0: '#FCD988'}, only_keys=None,
    max_position_to_plot = 200, add_to_title = '', show_num_actives = False):
        df_predictions = self.df_predictions

        # If only a subset of keys is requested
        if type(only_keys) is list and len(only_keys) > 0:
            keys_exist = np.all([ i in df_predictions.keys() for i in only_keys])
            assert keys_exist, 'Requested keys were not found.'
            df_predictions = {key: df_predictions[key] for key in only_keys}


        for key, y_pred in df_predictions.iteritems():
            order = np.argsort(-y_pred)
            y_pred_ord = y_pred[order]
            y_true = self.y_true[order]
            title = F'\n{add_to_title} {key} \n'
            if self.N > max_position_to_plot:
                y_true = y_true[:max_position_to_plot]
                # Get number of actives up to that position
                n_actives = (y_true[y_true == 1]).sum()
                title =  F'\n{add_to_title} {key} [first {max_position_to_plot} positions, {n_actives}/{self.n} actives found]\n'
            colors_array = [colors[i] for i in y_true]

            sns.palplot(sns.color_palette(colors_array))
            plt.title(title, fontsize=100)
            plt.show()

    # Formating metrics
    def format_metric_results(self, metric_name='roc_auc', only_keys=None,
                              rounded = 3, transposed = True, as_dataframe = True, **kwargs):
        if metric_name not in self.available_metrics:
            raise ValueError(F'Metric {metric_name} is not available. ' + 
                  F'Available metrics are:\n{self.available_metrics.keys()}')
        metric = self.available_metrics[metric_name]
        dic_results = {}

        df_predictions = self.df_predictions

        # If only a subset of keys is requested
        if type(only_keys) is list and len(only_keys) > 0:
            keys_exist = np.all([ i in df_predictions.keys() for i in only_keys])
            assert keys_exist, 'Requested keys were not found.'
            df_predictions = {key: df_predictions[key] for key in only_keys}


        for key, y_pred in df_predictions.iteritems():
            dic_results[key] = metric(y_pred, **kwargs)
        if as_dataframe:
            df = pd.DataFrame(dic_results, index = [metric_name.upper().replace('_', ' ')])
            df = df.T if transposed else df
            return df.round(rounded)
        else:
            return dic_results

#def plot_grid(list_of_plots, metric_to_plot, )