%include "sseutils.nasm"

section .data
section .bss
section .text

global sse

fob2        equ 16
dim_source  equ 12
distanze    equ 8



sse:
        push ebp
		mov	ebp, esp
		pushad   ;inizio
       
        mov     eax,  [ebp+distanze]        ; eax <- puntatore al vettore v delle distanze
        mov     ebx,  [ebp+dim_source]      ; ebx <- dimensione del vettore
        mov     ecx,  [ebp+fob2]            ; risultato(scalare) da restituire pari al quadrato degli elementi di 
                                            ; v sommati tra loro

        xor     esi,  esi                   ; i=0
        xorps   xmm4, xmm4                  ; xmm4=0

cicloQ: movaps  xmm0, [eax+esi*4]           ; p=4
        movaps  xmm1, [eax+esi*4+16]
        movaps  xmm2, [eax+esi*4+32]
        movaps  xmm3, [eax+esi*4+48]
        mulps   xmm0, xmm0
        mulps   xmm1, xmm1
        mulps   xmm2, xmm2
        mulps   xmm3, xmm3
        addps   xmm4, xmm1
        addps   xmm4, xmm2
        addps   xmm4, xmm3
        addps   xmm4, xmm0
        add     esi,  16
        cmp     esi,  ebx                  ; se esi<ebx -> i<dim_source
        jl cicloQ   

; todo - cicloR: 

        haddps xmm4, xmm4                  ;lezione 8
        haddps xmm4, xmm4
        movss [ecx], xmm4  





fine:   popad ;fine
		mov	esp, ebp
		pop	ebp
		ret
		
        
