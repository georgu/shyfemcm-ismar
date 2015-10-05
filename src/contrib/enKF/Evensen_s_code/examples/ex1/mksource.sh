#!/bin/ksh
#touch m_empty.F90
ls modules.F90 mod_*.F90 m_*.F90 |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR>0) {
       printf("MODULES = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'

#touch m_empty.F
ls mod_*.F m_*.F |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR>0) {
       printf("MODULES77 = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'

ls *.H |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR > 0) {
       printf("INC1 = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'


ls *.F90 | sed -e '/^modules/d' -e '/^mod_/d' -e '/^m_/d' -e '/^p_/d' |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR>0) {
       printf("F90FILES = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'


ls  *.F | sed -e '/^modules/d' -e '/^mod_/d' -e '/^m_/d' -e '/^p_/d' |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR>0) {
       printf("F77FILES = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'

