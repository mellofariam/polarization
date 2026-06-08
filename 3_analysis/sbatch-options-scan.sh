for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_contacts-scan.sh
    done
done 


for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_distance-histogram-scan.sh
    done
done 


for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_focused-pairwise-distance-distribution-scan.sh
    done
done


for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_energies-scan.sh
    done
done

for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_lamina-contact-scan.sh
    done
done


for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_msd-scan.sh
    done
done


for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_polarization-scan.sh
    done
done


for strength in 1 5 8; do 
    for condition in complete nucleolus; do
        sbatch --export=condition=$condition,strength=$strength sub_sasa-scan.sh
    done
done

