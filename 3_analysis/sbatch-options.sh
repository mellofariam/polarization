for condition in control complete isolation lamina nucleolus nuclear-bodies; do
    sbatch --export=chr_number=one,condition=$condition,chr=10 sub_msd.sh
done

for chr in 10 11; do
    for condition in control complete isolation lamina nucleolus nuclear-bodies; do
        sbatch --export=chr_number=two,condition=$condition,chr=$chr sub_msd.sh
    done
done

for condition in complete control lamina nuclear-bodies; do
    sbatch --export=chr_number=one,condition=$condition,chr=10 sub_msd.sh
done

for condition in nucleolus; do
    sbatch --export=chr_number=one,condition=$condition,chr=10 sub_msd.sh
done

for condition in control complete lamina nucleolus nuclear-bodies; do
    sbatch --export=chr_number=one,condition=$condition,chr=10 sub_msd.sh
done