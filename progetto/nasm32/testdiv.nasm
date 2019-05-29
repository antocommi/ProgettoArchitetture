extern printf

section .data
d db '%d '

section .bss
section .text

global testdiv

testdiv:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi	;inizio

        mov eax, 2
        mov ebx, 4
        mov ecx, 6
        mov edx, 8
        cdq

        pushad
        push eax
        push d
        call printf
        add esp, 8
        popad

        pushad
        push ebx
        push d
        call printf
        add esp, 8
        popad

        pushad
        push ecx
        push d
        call printf
        add esp, 8
        popad

        pushad
        push edx
        push d
        call printf
        add esp, 8
        popad

        idiv dword [ebp+8]

        pushad
        push eax
        push d
        call printf
        add esp, 8
        popad

        pushad
        push ebx
        push d
        call printf
        add esp, 8
        popad

        pushad
        push ecx
        push d
        call printf
        add esp, 8
        popad

        pushad
        push edx
        push d
        call printf
        add esp, 8
        popad

        pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret