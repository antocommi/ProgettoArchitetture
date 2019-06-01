%include "sseutils.nasm"

section .data
section .bss
section .text

global sse

fob2 equ 16
dim_source equ 12
distanza equ 8



sse:
        push ebp
		mov	ebp, esp
		pushad   ;inizio

       
        mov eax, [ebp+distanza]
        mov ebx, [ebp+dim_source]
        mov ecx, [ebp+fob2] ;risultato
        xor esi,esi
        xorps xmm4, xmm4
ciclo:  movaps xmm0, [eax+esi*4]
        movaps xmm1, [eax+esi*4+16]
        movaps xmm2, [eax+esi*4+32]
        movaps xmm3, [eax+esi*4+48]
        mulps xmm0, xmm0
        mulps xmm1, xmm1
        mulps xmm2, xmm2
        mulps xmm3, xmm3
        addps xmm4,xmm1
        addps xmm4,xmm2
        addps xmm4,xmm3
        addps xmm4,xmm0
        add esi,16
        cmp esi, ebx
        jl ciclo   

        haddps xmm4, xmm4   
        haddps xmm4, xmm4
        movss [ecx], xmm4  





fine:   popad ;fine
		mov	esp, ebp
		pop	ebp
		ret
		
        
