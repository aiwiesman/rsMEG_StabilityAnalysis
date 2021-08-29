function [median_ICCs,excellent_stability_times,good_stability_times,moderate_stability_times] = compute_temporalstability(data_array,num_perms,epoch_length)

%PURPOSE:           Compute temporal stability per each combination of region & frequency (ICC as a function of data length) as in: Wiesman AI, da Silva Castanheira J, Baillet S. (2021). Stability of spectral estimates in resting-state magnetoencephalography: recommendations for minimal data duration with neuroanatomical specificity.
%
%REQUIRED INPUTS:   data_array: 4-D array of data to be tested, with dimensions [region x frequency x participant x epoch]
%                   num_perms: the number of random permutations of the epoch order
%                   epoch_length: in seconds, for calculation of stability times
%
%OUTPUTS:           median_ICCs: 3-D array of ICCs (median over num_perms), with dimensions [region x frequency x epoch]
%                   excellent/good/moderate_stability_times: 2-D arrays representing the first time sample (in seconds) at which each [region x frequency] exceeded widely accepted ICC cutoffs
%
%NOTES:             dependency:  Arash Salarian (2021). Intraclass Correlation Coefficient (ICC) (https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc), MATLAB Central File Exchange. Retrieved August 17, 2021.
%                   ICC cutoffs from: Koo, T. K., & Li, M. Y. (2016). A guideline of selecting and reporting intraclass correlation coefficients for reliability research. Journal of chiropractic medicine, 15(2), 155-163.
%
%AUTHOR:            Alex I. Wiesman, neuroSPEED lab, Montreal Neurological Institute
%VERSION HISTORY:   08/05/2021  v1: First working version of program
%
%LICENSE:           This software is distributed under the terms of the GNU General Public License as published by the Free Software Foundation. Further details on the GPLv3 license can be found at http://www.gnu.org/copyleft/gpl.html.
%                   FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS DO NOT MAKE ANY WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.

for k = 1:num_perms %for each permutation, randomly shuffle data_array epoch order (dimension 4) then compute ICCs over data averaged across all possible independent epoch numbers from 1 to nepochs/2
    
    rand_dat = data_array(:,:,:,randperm(size(data_array,4)));  %shuffle epoch order
    fprintf('Running randomization %d/%d...\n',k,num_perms);
    
    for split = 1:size(rand_dat,4)/2    %%for each split from 1 to nepochs/2, average over epochs and compute ICC per each region and frequency
        first_section = mean(rand_dat(:,:,:,1:split),4);
        second_section = mean(rand_dat(:,:,:,(end-split)+1:end),4);
        for i = 1:size(rand_dat,1)
            for ii = 1:size(rand_dat,2)
                [r(i,ii,split,k), ~, ~, ~, ~, ~, ~] = ICC(squeeze([first_section(i,ii,:),second_section(i,ii,:)])', 'A-1', .05, .8);
            end
        end
    end
end

median_ICCs = median(r,4);  %%calculate median of ICCs across all randomizations

for i = 1:size(median_ICCs,1)   %%find first split-indices where median ICCs > .5 (moderate reliability)
    for ii = 1:size(median_ICCs,2)
        if any(median_ICCs(i,ii,:) > .50)
            moderate_split_ind_median(i,ii) = find(median_ICCs(i,ii,:) > .50,1,'first');
        else
            moderate_split_ind_median(i,ii) = NaN;
        end
    end
end
moderate_stability_times = moderate_split_ind_median.*epoch_length;

for i = 1:size(median_ICCs,1)   %%find first split-indices where median ICCs > .75 (good reliability)
    for ii = 1:size(median_ICCs,2)
        if any(median_ICCs(i,ii,:) > .75)
            good_split_ind_median(i,ii) = find(median_ICCs(i,ii,:) > .75,1,'first');
        else
            good_split_ind_median(i,ii) = NaN;
        end
    end
end
good_stability_times = good_split_ind_median.*epoch_length;

for i = 1:size(median_ICCs,1)   %%find first split-indices where median ICCs > .90 (excellent reliability)
    for ii = 1:size(median_ICCs,2)
        if any(median_ICCs(i,ii,:) > .90)
            excellent_split_ind_median(i,ii) = find(median_ICCs(i,ii,:) > .90,1,'first');
        else
            excellent_split_ind_median(i,ii) = NaN;
        end
    end
end
excellent_stability_times = excellent_split_ind_median.*epoch_length;
end