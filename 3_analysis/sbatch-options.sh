for condition in control complete lamina nucleolus; do
    sbatch --export=condition=$condition sub_msd.sh
done

for condition in control complete lamina nucleolus; do
    sbatch --export=condition=$condition sub_polarization.sh
done

for condition in complete nucleolus lamina control; do
    sbatch --export=condition=$condition sub_distance-histogram.sh
done

for condition in complete nucleolus lamina control; do
    sbatch --export=condition=$condition sub_focused-pairwise-distance-distribution.sh
done